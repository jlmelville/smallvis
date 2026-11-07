#include "bh.h"
#include <Rcpp.h>

using namespace Rcpp;

namespace {

struct CauchyKernel {
  double weight(double d2) const { return 1.0 / (1.0 + d2); }
  double slope(double d2) const { return weight(d2); }
  double log_weight(double d2) const { return -std::log(1.0 + d2); }
};

struct HeavyTailKernel {
  explicit HeavyTailKernel(double alpha)
      : alpha(alpha), exponent(-1.0 / alpha) {}

  double weight(double d2) const { return std::exp(log_weight(d2)); }
  double slope(double d2) const { return 1.0 / (1.0 + alpha * d2); }
  double log_weight(double d2) const {
    return exponent * std::log1p(alpha * d2);
  }

  double alpha;
  double exponent;
};

template <typename Accumulate>
void visit_bh(const smallvis::QuadTree::Node &node, double x, double y,
              double theta2, bool contains_query, Accumulate &accumulate) {
  if (node.num_points == 0) {
    return;
  }
  if (node.terminal_points) {
    for (const auto &point : *node.terminal_points) {
      const bool is_query = contains_query && point.x == x && point.y == y;
      const double dx = point.x - x;
      const double dy = point.y - y;
      accumulate(point.count - (is_query ? 1 : 0), dx, dy, dx * dx + dy * dy);
    }
    return;
  }

  const double dx = node.center_of_mass_x - x;
  const double dy = node.center_of_mass_y - y;
  const double d2 = dx * dx + dy * dy;
  if (node.is_summary(d2, theta2, contains_query)) {
    accumulate(node.num_points - (contains_query ? 1 : 0), dx, dy, d2);
    return;
  }

  const auto query_quadrant = node.find_quadrant(x, y);
  for (std::size_t qid = 0; qid < 4; ++qid) {
    visit_bh(*node.children[qid], x, y, theta2,
             contains_query && qid == query_quadrant, accumulate);
  }
}

template <bool Normalize, typename Kernel>
NumericMatrix
bh_gradient(const IntegerVector &indices, const IntegerVector &indptr,
            const NumericVector &P_data, const NumericMatrix &embedding,
            double theta, double eps, std::size_t n_threads,
            const Kernel &kernel, double repulsion) {
  if (embedding.ncol() != 2) {
    stop("Barnes-Hut requires two output dimensions");
  }

  const auto rows = as<std::vector<std::size_t>>(indices);
  const auto offsets = as<std::vector<std::size_t>>(indptr);
  const auto probabilities = as<std::vector<double>>(P_data);
  const auto Y = as<std::vector<double>>(transpose(embedding));
  const std::size_t n = Y.size() / 2;
  const smallvis::QuadTree tree(Y);
  const double theta2 = theta * theta;
  std::vector<double> gradient(Y.size(), 0.0);
  std::vector<double> sums(smallvis::worker_count(n, n_threads), 0.0);

  smallvis::parallel_for(
      n, n_threads,
      [&](std::size_t start, std::size_t end, std::size_t thread_id) {
        double sum = 0.0;
        for (std::size_t i = start; i < end; ++i) {
          double &gx = gradient[2 * i];
          double &gy = gradient[2 * i + 1];
          auto accumulate = [&](std::size_t count, double dx, double dy,
                                double d2) {
            const double mass = count * kernel.weight(d2);
            if (Normalize) {
              sum += mass;
            }
            const double force = mass * kernel.slope(d2);
            gx += force * dx;
            gy += force * dy;
          };
          visit_bh(*tree.root, Y[2 * i], Y[2 * i + 1], theta2, true,
                   accumulate);
        }
        if (Normalize) {
          sums[thread_id] = sum;
        }
      });

  if (Normalize) {
    const double Z = std::accumulate(sums.begin(), sums.end(), eps);
    for (auto &value : gradient) {
      value /= Z;
    }
  } else {
    for (auto &value : gradient) {
      value *= repulsion;
    }
  }

  smallvis::parallel_for(
      n, n_threads, [&](std::size_t start, std::size_t end, std::size_t) {
        for (std::size_t i = start; i < end; ++i) {
          for (std::size_t k = offsets[i]; k < offsets[i + 1]; ++k) {
            const std::size_t j = rows[k];
            const double dx = Y[2 * i] - Y[2 * j];
            const double dy = Y[2 * i + 1] - Y[2 * j + 1];
            const double force =
                probabilities[k] * kernel.slope(dx * dx + dy * dy);
            gradient[2 * i] += force * dx;
            gradient[2 * i + 1] += force * dy;
          }
        }
      });

  NumericMatrix result(2, embedding.nrow(), gradient.begin());
  return transpose(result);
}

// Sum P * log(w / Z) for normalized kernels, or P * log(w) - lambda * Z
// for elastic embedding. Both sums exclude self interactions.
template <typename Kernel>
double bh_score(const IntegerVector &indices, const IntegerVector &indptr,
                const NumericVector &P_data, const NumericMatrix &embedding,
                double theta, double eps, std::size_t n_threads,
                const Kernel &kernel, bool normalize, double repulsion) {
  if (embedding.ncol() != 2) {
    stop("Barnes-Hut requires two output dimensions");
  }

  const auto rows = as<std::vector<std::size_t>>(indices);
  const auto offsets = as<std::vector<std::size_t>>(indptr);
  const auto probabilities = as<std::vector<double>>(P_data);
  const auto Y = as<std::vector<double>>(transpose(embedding));
  const std::size_t n = Y.size() / 2;
  const smallvis::QuadTree tree(Y);
  const double theta2 = theta * theta;
  std::vector<double> sums(smallvis::worker_count(n, n_threads), 0.0);

  smallvis::parallel_for(
      n, n_threads,
      [&](std::size_t start, std::size_t end, std::size_t thread_id) {
        double sum = 0.0;
        for (std::size_t i = start; i < end; ++i) {
          auto accumulate = [&](std::size_t count, double, double, double d2) {
            sum += count * kernel.weight(d2);
          };
          visit_bh(*tree.root, Y[2 * i], Y[2 * i + 1], theta2, true,
                   accumulate);
        }
        sums[thread_id] = sum;
      });

  const double Z =
      std::accumulate(sums.begin(), sums.end(), normalize ? eps : 0.0);
  const double log_offset = normalize ? -std::log(Z) : 0.0;
  std::fill(sums.begin(), sums.end(), 0.0);
  smallvis::parallel_for(
      n, n_threads,
      [&](std::size_t start, std::size_t end, std::size_t thread_id) {
        double score = 0.0;
        for (std::size_t i = start; i < end; ++i) {
          for (std::size_t k = offsets[i]; k < offsets[i + 1]; ++k) {
            const std::size_t j = rows[k];
            const double dx = Y[2 * i] - Y[2 * j];
            const double dy = Y[2 * i + 1] - Y[2 * j + 1];
            score += probabilities[k] *
                     (log_offset + kernel.log_weight(dx * dx + dy * dy));
          }
        }
        sums[thread_id] = score;
      });
  return std::accumulate(sums.begin(), sums.end(),
                         normalize ? 0.0 : -repulsion * Z);
}

} // namespace

// [[Rcpp::export]]
NumericMatrix bh_tsne_gradient_cpp(const IntegerVector &indices,
                                   const IntegerVector &indptr,
                                   const NumericVector &P_data,
                                   const NumericMatrix &embedding,
                                   double theta = 0.5, double eps = 1e-16,
                                   std::size_t n_threads = 1) {
  return bh_gradient<true>(indices, indptr, P_data, embedding, theta, eps,
                           n_threads, CauchyKernel{}, 0.0);
}

// [[Rcpp::export]]
double bh_plogq_cpp(const IntegerVector &indices, const IntegerVector &indptr,
                    const NumericVector &P_data, const NumericMatrix &embedding,
                    double theta = 0.5, double eps = 1e-16,
                    std::size_t n_threads = 1) {
  return bh_score(indices, indptr, P_data, embedding, theta, eps, n_threads,
                  CauchyKernel{}, true, 0.0);
}

// [[Rcpp::export]]
NumericMatrix bh_sne_gradient_cpp(const IntegerVector &indices,
                                  const IntegerVector &indptr,
                                  const NumericVector &P_data,
                                  const NumericMatrix &embedding, double theta,
                                  double eps, int n_threads, double alpha,
                                  bool normalize, double lambda) {
  if (normalize && alpha != 1.0) {
    return bh_gradient<true>(indices, indptr, P_data, embedding, theta, eps,
                             n_threads, HeavyTailKernel(alpha), 0.0);
  }
  if (normalize) {
    return bh_tsne_gradient_cpp(indices, indptr, P_data, embedding, theta, eps,
                                n_threads);
  }
  return bh_gradient<false>(indices, indptr, P_data, embedding, theta, eps,
                            n_threads, CauchyKernel{}, lambda);
}

// [[Rcpp::export]]
double bh_sne_score_cpp(const IntegerVector &indices,
                        const IntegerVector &indptr,
                        const NumericVector &P_data,
                        const NumericMatrix &embedding, double theta,
                        double eps, int n_threads, double alpha, bool normalize,
                        double lambda) {
  if (normalize && alpha != 1.0) {
    return bh_score(indices, indptr, P_data, embedding, theta, eps, n_threads,
                    HeavyTailKernel(alpha), true, 0.0);
  }
  return bh_score(indices, indptr, P_data, embedding, theta, eps, n_threads,
                  CauchyKernel{}, normalize, lambda);
}
