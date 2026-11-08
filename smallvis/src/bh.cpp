#include "bh.h"
#include <Rcpp.h>

using namespace Rcpp;

namespace {

struct CauchyKernel {
  double weight(double d2) const { return 1.0 / (1.0 + d2); }
  double slope(double d2) const { return weight(d2); }
  double log_weight(double d2) const { return -std::log(1.0 + d2); }
  double repulsive_force(double, double w, double mass) const {
    return mass * w;
  }
  double score_mass(std::size_t count, double d2) const {
    return count * weight(d2);
  }
};

struct HeavyTailKernel {
  explicit HeavyTailKernel(double alpha)
      : alpha(alpha), exponent(-1.0 / alpha) {}

  double weight(double d2) const { return std::exp(log_weight(d2)); }
  double slope(double d2) const { return 1.0 / (1.0 + alpha * d2); }
  double log_weight(double d2) const {
    return exponent * std::log1p(alpha * d2);
  }
  double repulsive_force(double d2, double, double mass) const {
    return mass * slope(d2);
  }
  double score_mass(std::size_t count, double d2) const {
    return count * weight(d2);
  }

  double alpha;
  double exponent;
};

struct LargeVisKernel : CauchyKernel {
  LargeVisKernel(double gr_eps, double eps) : gr_eps(gr_eps), eps(eps) {}

  double repulsive_force(double d2, double, double mass) const {
    if (mass == 0) {
      return 0; // Self-only leaf.
    }
    if (d2 == 0 && gr_eps == 0) {
      throw std::domain_error(
          "LargeVis with gr_eps = 0 has undefined repulsion at zero distance");
    }
    // w^2 / (1 + (gr_eps - 1) * w) = w / (d2 + gr_eps).
    return mass / (d2 + gr_eps);
  }

  double score_mass(std::size_t count, double d2) const {
    return count == 0 ? 0.0 : -(count * std::log1p(eps - weight(d2)));
  }

  double gr_eps;
  double eps;
};

template <bool Student> struct UmapKernel {
  UmapKernel(double a, double b, double gr_eps, double eps)
      : a(a), b(b), stabilization(gr_eps + eps), eps(eps) {}

  double weight(double d2) const {
    return Student ? CauchyKernel{}.weight(d2)
                   : 1.0 / (1.0 + a * std::pow(d2, b));
  }

  double slope(double d2) const {
    // For t-UMAP, w / (1 - w + eps + gr_eps), with w = 1 / (1 + d2).
    return Student ? 1.0 / (d2 + stabilization * (1.0 + d2))
                   : b / (d2 + stabilization);
  }
  double repulsive_force(double d2, double, double mass) const {
    return mass * slope(d2);
  }
  double log_weight(double d2) const {
    const double w = weight(d2);
    return std::log(w + eps) - std::log1p(eps - w);
  }
  double score_mass(std::size_t count, double d2) const {
    // Subtract the cost at w = 0, so the diagonal contributes zero.
    return count == 0
               ? 0.0
               : count * (std::log1p(eps) - std::log1p(eps - weight(d2)));
  }

  double a;
  double b;
  double stabilization;
  double eps;
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
double bh_repulsion(const smallvis::QuadTree &tree,
                    const std::vector<double> &Y, double theta, double eps,
                    std::size_t n_threads, const Kernel &kernel,
                    std::vector<double> &gradient) {
  const std::size_t n = Y.size() / 2;
  const double theta2 = theta * theta;
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
            const double w = kernel.weight(d2);
            const double mass = count * w;
            if (Normalize) {
              sum += mass;
            }
            const double force = kernel.repulsive_force(d2, w, mass);
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

  return Normalize ? std::accumulate(sums.begin(), sums.end(), eps) : 0.0;
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
  std::vector<double> gradient(Y.size(), 0.0);
  const double Z =
      bh_repulsion<Normalize>(tree, Y, theta, eps, n_threads, kernel, gradient);

  if (Normalize) {
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

// The kernel supplies either normalization mass or repulsive cost.
template <typename Kernel>
double bh_score_mass(const std::vector<double> &Y, double theta, double eps,
                     std::size_t n_threads, const Kernel &kernel,
                     std::vector<double> &sums) {
  const std::size_t n = Y.size() / 2;
  const smallvis::QuadTree tree(Y);
  const double theta2 = theta * theta;

  smallvis::parallel_for(
      n, n_threads,
      [&](std::size_t start, std::size_t end, std::size_t thread_id) {
        double sum = 0.0;
        for (std::size_t i = start; i < end; ++i) {
          auto accumulate = [&](std::size_t count, double, double, double d2) {
            sum += kernel.score_mass(count, d2);
          };
          visit_bh(*tree.root, Y[2 * i], Y[2 * i + 1], theta2, true,
                   accumulate);
        }
        sums[thread_id] = sum;
      });

  return std::accumulate(sums.begin(), sums.end(), eps);
}

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
  std::vector<double> sums(smallvis::worker_count(n, n_threads), 0.0);
  const double Z =
      bh_score_mass(Y, theta, normalize ? eps : 0.0, n_threads, kernel, sums);
  const double log_offset = normalize ? -std::log(Z) : 0.0;
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

template <bool Gradient>
double jssne_edge_sum(const IntegerVector &indices, const IntegerVector &indptr,
                      const NumericVector &P_data, const std::vector<double> &Y,
                      double Z, std::size_t n_threads,
                      std::vector<double> *gradient) {
  const auto rows = as<std::vector<std::size_t>>(indices);
  const auto offsets = as<std::vector<std::size_t>>(indptr);
  const auto probabilities = as<std::vector<double>>(P_data);
  const std::size_t n = Y.size() / 2;
  const CauchyKernel kernel;
  std::vector<double> sums(smallvis::worker_count(n, n_threads), 0.0);
  smallvis::parallel_for(
      n, n_threads,
      [&](std::size_t start, std::size_t end, std::size_t thread_id) {
        double sum = 0.0;
        for (std::size_t i = start; i < end; ++i) {
          for (std::size_t k = offsets[i]; k < offsets[i + 1]; ++k) {
            const std::size_t j = rows[k];
            const double dx = Y[2 * i] - Y[2 * j];
            const double dy = Y[2 * i + 1] - Y[2 * j + 1];
            const double w = kernel.weight(dx * dx + dy * dy);
            const double q = w / Z;
            const double p = probabilities[k];
            // log(2Q / (P + Q)) = log(2) - log1p(P / Q).
            // The correction vanishes off the input graph.
            const double correction = std::log1p(p / q);
            sum += q * correction;
            if (Gradient) {
              const double force = w * w * correction;
              (*gradient)[2 * i] += force * dx;
              (*gradient)[2 * i + 1] += force * dy;
            } else {
              sum += p * std::log(0.5 * (p + q));
            }
          }
        }
        sums[thread_id] = sum;
      });
  return std::accumulate(sums.begin(), sums.end(), 0.0);
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

// [[Rcpp::export]]
NumericMatrix bh_largevis_gradient_cpp(const IntegerVector &indices,
                                       const IntegerVector &indptr,
                                       const NumericVector &P_data,
                                       const NumericMatrix &embedding,
                                       double theta, int n_threads,
                                       double gamma, double gr_eps) {
  return bh_gradient<false>(indices, indptr, P_data, embedding, theta, 0.0,
                            n_threads, LargeVisKernel(gr_eps, 0.0), gamma);
}

// [[Rcpp::export]]
double bh_largevis_score_cpp(const IntegerVector &indices,
                             const IntegerVector &indptr,
                             const NumericVector &P_data,
                             const NumericMatrix &embedding, double theta,
                             int n_threads, double gamma, double eps) {
  const double score =
      bh_score(indices, indptr, P_data, embedding, theta, 0.0, n_threads,
               LargeVisKernel(0.0, eps), false, gamma);
  // The dense objective includes the constant diagonal, where w = 0.
  return score + embedding.nrow() * gamma * std::log1p(eps);
}

// [[Rcpp::export]]
NumericMatrix bh_tumap_gradient_cpp(const IntegerVector &indices,
                                    const IntegerVector &indptr,
                                    const NumericVector &P_data,
                                    const NumericMatrix &embedding,
                                    double theta, int n_threads, double gr_eps,
                                    double eps) {
  return bh_gradient<false>(indices, indptr, P_data, embedding, theta, 0.0,
                            n_threads, UmapKernel<true>(1.0, 1.0, gr_eps, eps),
                            1.0);
}

// [[Rcpp::export]]
double bh_tumap_score_cpp(const IntegerVector &indices,
                          const IntegerVector &indptr,
                          const NumericVector &P_data,
                          const NumericMatrix &embedding, double theta,
                          int n_threads, double eps) {
  return bh_score(indices, indptr, P_data, embedding, theta, 0.0, n_threads,
                  UmapKernel<true>(1.0, 1.0, 0.0, eps), false, 1.0);
}

// [[Rcpp::export]]
NumericMatrix bh_umap_gradient_cpp(const IntegerVector &indices,
                                   const IntegerVector &indptr,
                                   const NumericVector &P_data,
                                   const NumericMatrix &embedding, double theta,
                                   int n_threads, double gr_eps, double eps,
                                   double a, double b) {
  return bh_gradient<false>(indices, indptr, P_data, embedding, theta, 0.0,
                            n_threads, UmapKernel<false>(a, b, gr_eps, eps),
                            1.0);
}

// [[Rcpp::export]]
double bh_umap_score_cpp(const IntegerVector &indices,
                         const IntegerVector &indptr,
                         const NumericVector &P_data,
                         const NumericMatrix &embedding, double theta,
                         int n_threads, double eps, double a, double b) {
  return bh_score(indices, indptr, P_data, embedding, theta, 0.0, n_threads,
                  UmapKernel<false>(a, b, 0.0, eps), false, 1.0);
}

// [[Rcpp::export]]
NumericMatrix bh_jssne_gradient_cpp(const IntegerVector &indices,
                                    const IntegerVector &indptr,
                                    const NumericVector &P_data,
                                    const NumericMatrix &embedding,
                                    double theta, int n_threads) {
  if (embedding.ncol() != 2) {
    stop("Barnes-Hut requires two output dimensions");
  }
  const auto Y = as<std::vector<double>>(transpose(embedding));
  const smallvis::QuadTree tree(Y);
  std::vector<double> repulsion(Y.size(), 0.0);
  const double Z = bh_repulsion<true>(tree, Y, theta, 0.0, n_threads,
                                      CauchyKernel{}, repulsion);
  std::vector<double> gradient(Y.size(), 0.0);
  const double edge_sum =
      jssne_edge_sum<true>(indices, indptr, P_data, Y, Z, n_threads, &gradient);
  const double scale = 2.0 / Z;
  for (std::size_t i = 0; i < gradient.size(); ++i) {
    gradient[i] = scale * (gradient[i] + edge_sum * repulsion[i]);
  }
  const NumericMatrix result(2, embedding.nrow(), gradient.begin());
  return transpose(result);
}

// [[Rcpp::export]]
double bh_jssne_score_cpp(const IntegerVector &indices,
                          const IntegerVector &indptr,
                          const NumericVector &P_data,
                          const NumericMatrix &embedding, double theta,
                          int n_threads) {
  if (embedding.ncol() != 2) {
    stop("Barnes-Hut requires two output dimensions");
  }
  const auto Y = as<std::vector<double>>(transpose(embedding));
  const std::size_t n = Y.size() / 2;
  std::vector<double> sums(smallvis::worker_count(n, n_threads), 0.0);
  const double Z =
      bh_score_mass(Y, theta, 0.0, n_threads, CauchyKernel{}, sums);
  return std::log(2.0) - jssne_edge_sum<false>(indices, indptr, P_data, Y, Z,
                                               n_threads, nullptr);
}
