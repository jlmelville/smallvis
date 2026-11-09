#include "bh.h"
#include <Rcpp.h>
#include <type_traits>

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

struct GaussianKernel {
  double weight(double d2) const {
    // Only self-only leaves can have d2 below the shift.
    return std::exp(std::min(0.0, log_weight(d2)));
  }
  double slope(double) const { return 1.0; }
  double log_weight(double d2) const { return shift - d2; }
  double repulsive_force(double, double, double mass) const { return mass; }
  double score_mass(std::size_t count, double d2) const {
    return count * weight(d2);
  }

  double shift = 0.0;
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

template <typename Kernel>
bool rescale_weights(Kernel &, const smallvis::QuadTree &,
                     const std::vector<double> &, double, double) {
  return false;
}

bool rescale_weights(GaussianKernel &kernel, const smallvis::QuadTree &tree,
                     const std::vector<double> &Y, double theta, double Z) {
  // Shift only when the unscaled sum is close to underflow. Use the same
  // representatives as the normalization and exclude self-only leaves.
  const double min_mass = std::numeric_limits<double>::min() /
                          std::numeric_limits<double>::epsilon();
  if (Z >= min_mass) {
    return false;
  }
  double min_d2 = std::numeric_limits<double>::infinity();
  auto accumulate = [&](std::size_t count, double, double, double d2) {
    if (count != 0) {
      min_d2 = std::min(min_d2, d2);
    }
  };
  for (std::size_t i = 0; i < Y.size(); i += 2) {
    visit_bh(*tree.root, Y[i], Y[i + 1], theta * theta, true, accumulate);
  }
  if (!std::isfinite(min_d2)) {
    throw std::domain_error("Barnes-Hut SSNE requires a finite pair distance");
  }
  kernel.shift = min_d2;
  return true;
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

// Accumulate the sparse attractive force into an interleaved gradient.
template <typename Index, typename Kernel>
void sparse_attraction(const Index *rows, const Index *offsets,
                       const double *probabilities,
                       const std::vector<double> &Y, std::size_t n_threads,
                       const Kernel &kernel, std::vector<double> &gradient) {
  const std::size_t n = Y.size() / 2;
  smallvis::parallel_for(
      n, n_threads, [&](std::size_t start, std::size_t end, std::size_t) {
        for (std::size_t i = start; i < end; ++i) {
          for (Index k = offsets[i]; k < offsets[i + 1]; ++k) {
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
}

template <bool Normalize, typename Kernel>
NumericMatrix
bh_gradient(const IntegerVector &indices, const IntegerVector &indptr,
            const NumericVector &P_data, const NumericMatrix &embedding,
            double theta, double eps, std::size_t n_threads, Kernel kernel,
            double repulsion) {
  if (embedding.ncol() != 2) {
    stop("Barnes-Hut requires two output dimensions");
  }

  const auto rows = as<std::vector<std::size_t>>(indices);
  const auto offsets = as<std::vector<std::size_t>>(indptr);
  const auto probabilities = as<std::vector<double>>(P_data);
  const auto Y = as<std::vector<double>>(transpose(embedding));
  const smallvis::QuadTree tree(Y);
  std::vector<double> gradient(Y.size(), 0.0);
  double Z =
      bh_repulsion<Normalize>(tree, Y, theta, eps, n_threads, kernel, gradient);

  if (Normalize) {
    if (rescale_weights(kernel, tree, Y, theta, Z)) {
      std::fill(gradient.begin(), gradient.end(), 0.0);
      Z = bh_repulsion<Normalize>(tree, Y, theta, eps, n_threads, kernel,
                                  gradient);
    }
    for (auto &value : gradient) {
      value /= Z;
    }
  } else {
    for (auto &value : gradient) {
      value *= repulsion;
    }
  }

  sparse_attraction(rows.data(), offsets.data(), probabilities.data(), Y,
                    n_threads, kernel, gradient);

  NumericMatrix result(2, embedding.nrow(), gradient.begin());
  return transpose(result);
}

// The kernel supplies either normalization mass or repulsive cost.
template <typename Kernel>
double bh_score_mass(const smallvis::QuadTree &tree,
                     const std::vector<double> &Y, double theta, double eps,
                     std::size_t n_threads, const Kernel &kernel,
                     std::vector<double> &sums) {
  const std::size_t n = Y.size() / 2;
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

template <typename Index, typename Kernel>
double sparse_log_weight_sum(const Index *rows, const Index *offsets,
                             const double *probabilities,
                             const std::vector<double> &Y,
                             std::size_t n_threads, const Kernel &kernel,
                             double log_offset, double initial) {
  const std::size_t n = Y.size() / 2;
  std::vector<double> sums(smallvis::worker_count(n, n_threads), 0.0);
  smallvis::parallel_for(
      n, n_threads,
      [&](std::size_t start, std::size_t end, std::size_t thread_id) {
        double score = 0.0;
        for (std::size_t i = start; i < end; ++i) {
          for (Index k = offsets[i]; k < offsets[i + 1]; ++k) {
            const std::size_t j = rows[k];
            const double dx = Y[2 * i] - Y[2 * j];
            const double dy = Y[2 * i + 1] - Y[2 * j + 1];
            score += probabilities[k] *
                     (log_offset + kernel.log_weight(dx * dx + dy * dy));
          }
        }
        sums[thread_id] = score;
      });
  return std::accumulate(sums.begin(), sums.end(), initial);
}

template <typename Kernel>
double bh_score(const IntegerVector &indices, const IntegerVector &indptr,
                const NumericVector &P_data, const NumericMatrix &embedding,
                double theta, double eps, std::size_t n_threads, Kernel kernel,
                bool normalize, double repulsion) {
  if (embedding.ncol() != 2) {
    stop("Barnes-Hut requires two output dimensions");
  }

  const auto rows = as<std::vector<std::size_t>>(indices);
  const auto offsets = as<std::vector<std::size_t>>(indptr);
  const auto probabilities = as<std::vector<double>>(P_data);
  const auto Y = as<std::vector<double>>(transpose(embedding));
  const std::size_t n = Y.size() / 2;
  const smallvis::QuadTree tree(Y);
  std::vector<double> sums(smallvis::worker_count(n, n_threads), 0.0);
  const double offset = normalize ? eps : 0.0;
  double Z = bh_score_mass(tree, Y, theta, offset, n_threads, kernel, sums);
  if (normalize && rescale_weights(kernel, tree, Y, theta, Z)) {
    Z = bh_score_mass(tree, Y, theta, offset, n_threads, kernel, sums);
  }
  const double log_offset = normalize ? -std::log(Z) : 0.0;
  return sparse_log_weight_sum(rows.data(), offsets.data(),
                               probabilities.data(), Y, n_threads, kernel,
                               log_offset, normalize ? 0.0 : -repulsion * Z);
}

enum class Divergence { JensenShannon, Hellinger, ChiSquared };

template <Divergence Kind, bool Gradient>
double
divergence_edge_sum(const IntegerVector &indices, const IntegerVector &indptr,
                    const NumericVector &edge_data,
                    const std::vector<double> &Y, double Z,
                    std::size_t n_threads, std::vector<double> *gradient) {
  const auto rows = as<std::vector<std::size_t>>(indices);
  const auto offsets = as<std::vector<std::size_t>>(indptr);
  const auto edges = as<std::vector<double>>(edge_data);
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
            const double d2 = dx * dx + dy * dy;
            const double edge = edges[k];
            double force = 0.0;
            if (Kind == Divergence::ChiSquared) {
              // The input edge holds P^2, so this sums P^2 / w.
              sum += edge * (1.0 + d2);
              force = edge;
            } else {
              const double w = kernel.weight(d2);
              const double q = w / Z;
              if (Kind == Divergence::Hellinger) {
                // The input edge already holds sqrt(P).
                const double sqrt_q = std::sqrt(q);
                sum += edge * sqrt_q;
                if (Gradient) {
                  force = w * w * (edge / sqrt_q);
                }
              } else {
                // log(2Q / (P + Q)) = log(2) - log1p(P / Q).
                const double correction = std::log1p(edge / q);
                sum += q * correction;
                if (Gradient) {
                  force = w * w * correction;
                } else {
                  sum += edge * std::log(0.5 * (edge + q));
                }
              }
            }
            if (Gradient) {
              (*gradient)[2 * i] += force * dx;
              (*gradient)[2 * i + 1] += force * dy;
            }
          }
        }
        sums[thread_id] = sum;
      });
  return std::accumulate(sums.begin(), sums.end(), 0.0);
}

template <Divergence Kind>
NumericMatrix divergence_gradient_from_sums(
    const IntegerVector &indices, const IntegerVector &indptr,
    const NumericVector &edge_data, const std::vector<double> &Y, double Z,
    const std::vector<double> &repulsion, double repulsion_sign,
    std::size_t n_threads) {
  std::vector<double> gradient(Y.size(), 0.0);
  const double edge_sum = divergence_edge_sum<Kind, true>(
      indices, indptr, edge_data, Y, Z, n_threads, &gradient);
  const bool chi_squared = Kind == Divergence::ChiSquared;
  const double scale =
      (Kind == Divergence::JensenShannon ? 2.0 : 4.0) / (chi_squared ? 1.0 : Z);
  const double attraction_scale = chi_squared ? Z : 1.0;
  const double signed_edge_sum = repulsion_sign * edge_sum;
  for (std::size_t i = 0; i < gradient.size(); ++i) {
    gradient[i] = scale * (attraction_scale * gradient[i] +
                           signed_edge_sum * repulsion[i]);
  }
  const NumericMatrix result(2, Y.size() / 2, gradient.begin());
  return transpose(result);
}

template <Divergence Kind>
NumericMatrix bh_divergence_gradient(const IntegerVector &indices,
                                     const IntegerVector &indptr,
                                     const NumericVector &edge_data,
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
  return divergence_gradient_from_sums<Kind>(indices, indptr, edge_data, Y, Z,
                                             repulsion, 1.0, n_threads);
}

template <Divergence Kind>
double divergence_score_from_mass(const IntegerVector &indices,
                                  const IntegerVector &indptr,
                                  const NumericVector &edge_data,
                                  const std::vector<double> &Y, double Z,
                                  std::size_t n_threads) {
  const double edge_sum = divergence_edge_sum<Kind, false>(
      indices, indptr, edge_data, Y, Z, n_threads, nullptr);
  if (Kind == Divergence::JensenShannon) {
    return std::log(2.0) - edge_sum;
  }
  return Kind == Divergence::ChiSquared ? Z * edge_sum : edge_sum;
}

template <Divergence Kind>
double bh_divergence_score(const IntegerVector &indices,
                           const IntegerVector &indptr,
                           const NumericVector &edge_data,
                           const NumericMatrix &embedding, double theta,
                           int n_threads) {
  if (embedding.ncol() != 2) {
    stop("Barnes-Hut requires two output dimensions");
  }
  const auto Y = as<std::vector<double>>(transpose(embedding));
  const std::size_t n = Y.size() / 2;
  const smallvis::QuadTree tree(Y);
  std::vector<double> sums(smallvis::worker_count(n, n_threads), 0.0);
  const double Z =
      bh_score_mass(tree, Y, theta, 0.0, n_threads, CauchyKernel{}, sums);
  return divergence_score_from_mass<Kind>(indices, indptr, edge_data, Y, Z,
                                          n_threads);
}

struct KlRow {
  double mass = 0.0;
  double log_mass = 0.0;
  double forward_loss = 0.0;
  double force[2] = {};
  double log_force[2] = {};
};

template <bool Gradient, typename Kernel>
double bh_kl_tree_rows(const smallvis::QuadTree &tree,
                       const std::vector<double> &Y, double theta,
                       std::size_t n_threads, const Kernel &kernel,
                       std::vector<KlRow> &rows) {
  const double theta2 = theta * theta;
  smallvis::parallel_for(rows.size(), n_threads,
                         [&](std::size_t start, std::size_t end, std::size_t) {
                           for (std::size_t i = start; i < end; ++i) {
                             KlRow row;
                             auto accumulate = [&](std::size_t count, double dx,
                                                   double dy, double d2) {
                               const double mass = count * kernel.weight(d2);
                               if (mass == 0.0) {
                                 return;
                               }
                               const double log_w = kernel.log_weight(d2);
                               row.mass += mass;
                               row.log_mass += mass * log_w;
                               if (Gradient) {
                                 const double force = mass * kernel.slope(d2);
                                 const double direction[2] = {dx, dy};
                                 for (std::size_t axis = 0; axis < 2; ++axis) {
                                   const double component =
                                       force * direction[axis];
                                   row.force[axis] += component;
                                   row.log_force[axis] += component * log_w;
                                 }
                               }
                             };
                             visit_bh(*tree.root, Y[2 * i], Y[2 * i + 1],
                                      theta2, true, accumulate);
                             rows[i] = row;
                           }
                         });
  return std::accumulate(
      rows.begin(), rows.end(), 0.0,
      [](double sum, const KlRow &row) { return sum + row.mass; });
}

template <bool Gradient, typename Kernel>
SEXP mixed_kl_from_sums(const IntegerVector &indices,
                        const IntegerVector &indptr,
                        const NumericVector &log_ratio,
                        const NumericVector &p_excess,
                        const std::vector<double> &Y, const double *center,
                        std::vector<KlRow> &rows, double Z,
                        std::size_t n_threads, double eps, double lambda,
                        double p_mass, double plogp, const Kernel &kernel) {
  const auto neighbors = as<std::vector<std::size_t>>(indices);
  const auto offsets = as<std::vector<std::size_t>>(indptr);
  const auto logs = as<std::vector<double>>(log_ratio);
  const auto excess = as<std::vector<double>>(p_excess);
  const std::size_t n = Y.size() / 2;
  // Every off-diagonal pair starts at P = eps. Sparse edges correct log P
  // in the reverse term and P itself in the forward term.
  const double reverse_weight = (1.0 - lambda) / Z;
  std::vector<double> gradient(Gradient ? Y.size() : 0, 0.0);
  smallvis::parallel_for(
      n, n_threads, [&](std::size_t start, std::size_t end, std::size_t) {
        for (std::size_t i = start; i < end; ++i) {
          for (std::size_t k = offsets[i]; k < offsets[i + 1]; ++k) {
            const std::size_t j = neighbors[k];
            const double dx = Y[2 * i] - Y[2 * j];
            const double dy = Y[2 * i + 1] - Y[2 * j + 1];
            const double d2 = dx * dx + dy * dy;
            const double correction = kernel.weight(d2) * logs[k];
            rows[i].log_mass -= correction;
            if (Gradient) {
              const double force =
                  reverse_weight * correction * kernel.slope(d2) +
                  lambda * excess[k];
              gradient[2 * i] += force * dx;
              gradient[2 * i + 1] += force * dy;
            } else if (lambda != 0.0) {
              rows[i].forward_loss -= excess[k] * kernel.log_weight(d2);
            }
          }
        }
      });

  const double mean_log = std::accumulate(rows.begin(), rows.end(), 0.0,
                                          [](double sum, const KlRow &row) {
                                            return sum + row.log_mass;
                                          }) /
                          Z;
  const double forward_weight = lambda * p_mass / Z;
  const double floor_force = lambda * eps * n;
  double forward_loss = 0.0;
  double centered_d2 = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    const KlRow &row = rows[i];
    for (std::size_t axis = 0; axis < 2; ++axis) {
      const double displacement = Y[2 * i + axis] - center[axis];
      if (Gradient) {
        gradient[2 * i + axis] +=
            reverse_weight *
                (row.log_force[axis] - mean_log * row.force[axis]) +
            forward_weight * row.force[axis] + floor_force * displacement;
        gradient[2 * i + axis] *= 4.0;
      } else if (lambda != 0.0) {
        centered_d2 += displacement * displacement;
      }
    }
    forward_loss += row.forward_loss;
  }
  if (Gradient) {
    NumericMatrix result(2, n, gradient.begin());
    return transpose(result);
  }

  const double log_z = std::log(Z);
  const double reverse = mean_log - log_z - std::log(eps);
  if (lambda == 0.0) {
    return wrap(reverse);
  }
  // The Gaussian floor term uses the sum of all squared pair distances.
  const double floor_loss =
      eps * n * (2.0 * centered_d2 - (n - 1.0) * kernel.log_weight(0.0));
  const double forward = plogp + p_mass * log_z + forward_loss + floor_loss;
  return wrap(lambda * forward + (1.0 - lambda) * reverse);
}

template <bool Gradient, typename Kernel>
SEXP bh_mixed_kl(const IntegerVector &indices, const IntegerVector &indptr,
                 const NumericVector &log_ratio, const NumericVector &p_excess,
                 const NumericMatrix &embedding, double theta,
                 std::size_t n_threads, double eps, double lambda,
                 double p_mass, double plogp, Kernel kernel) {
  if (embedding.ncol() != 2) {
    stop("Barnes-Hut requires two output dimensions");
  }
  const auto Y = as<std::vector<double>>(transpose(embedding));
  const std::size_t n = Y.size() / 2;
  const smallvis::QuadTree tree(Y);
  std::vector<KlRow> rows(n);
  double Z = bh_kl_tree_rows<Gradient>(tree, Y, theta, n_threads, kernel, rows);
  if (rescale_weights(kernel, tree, Y, theta, Z)) {
    Z = bh_kl_tree_rows<Gradient>(tree, Y, theta, n_threads, kernel, rows);
  }
  if (!(Z > 0.0) || !std::isfinite(Z)) {
    stop("Barnes-Hut KL requires a finite positive normalizer");
  }

  const double center[2] = {tree.root->center_of_mass_x,
                            tree.root->center_of_mass_y};
  return mixed_kl_from_sums<Gradient>(indices, indptr, log_ratio, p_excess, Y,
                                      center, rows, Z, n_threads, eps, lambda,
                                      p_mass, plogp, kernel);
}

template <typename Evaluate>
auto with_bernoulli_kernel(const std::string &method, double eps, double gr_eps,
                           double a, double b, Evaluate evaluate) {
  if (method == "largevis") {
    return evaluate(LargeVisKernel(gr_eps, eps));
  }
  if (method == "tumap") {
    return evaluate(UmapKernel<true>(1.0, 1.0, gr_eps, eps));
  }
  if (method == "umap") {
    return evaluate(UmapKernel<false>(a, b, gr_eps, eps));
  }
  stop("Unknown FFT Bernoulli kernel");
}

// Smooth the radial kernel across the cutoff with a quintic transition.
// The local correction adds back exactly the part removed by this split.
template <typename Kernel> struct BernoulliGridKernel {
  BernoulliGridKernel(const Kernel &kernel, double radius, bool gradient)
      : kernel(kernel), cutoff2(radius * radius), gradient(gradient),
        at_cutoff(radius > 0.0 ? exact(cutoff2) : 0.0) {}

  double exact(double d2) const {
    if (gradient) {
      const double w = kernel.weight(d2);
      return kernel.repulsive_force(d2, w, w);
    }
    return kernel.score_mass(1, d2);
  }

  double near_fraction(double d2) const {
    const double t = d2 / cutoff2;
    const double complement = 1.0 - t;
    return complement * complement * complement * (1.0 + t * (3.0 + 6.0 * t));
  }

  double smooth(double d2) const {
    if (d2 >= cutoff2) {
      return exact(d2);
    }
    if (d2 == 0.0) {
      return at_cutoff;
    }
    const double fraction = near_fraction(d2);
    return at_cutoff * fraction + (1.0 - fraction) * exact(d2);
  }

  double correction(double d2) const {
    return near_fraction(d2) * (exact(d2) - at_cutoff);
  }

  const Kernel &kernel;
  double cutoff2;
  bool gradient;
  double at_cutoff;
};

// Bucket widths equal the correction radius, so only adjacent buckets can
// contain interacting pairs. Each worker updates only its own query points.
template <typename Kernel>
double bernoulli_near_pairs(const std::vector<double> &Y, double radius,
                            std::size_t n_threads,
                            const BernoulliGridKernel<Kernel> &split,
                            std::vector<double> &gradient) {
  if (radius == 0.0) {
    return 0.0;
  }
  const std::size_t n = Y.size() / 2;
  double xmin = Y[0], ymin = Y[1], xmax = xmin, ymax = ymin;
  for (std::size_t i = 1; i < n; ++i) {
    xmin = std::min(xmin, Y[2 * i]);
    xmax = std::max(xmax, Y[2 * i]);
    ymin = std::min(ymin, Y[2 * i + 1]);
    ymax = std::max(ymax, Y[2 * i + 1]);
  }
  const auto nx = static_cast<std::size_t>((xmax - xmin) / radius) + 1;
  const auto ny = static_cast<std::size_t>((ymax - ymin) / radius) + 1;
  std::vector<std::size_t> cell(n), offsets(nx * ny + 1, 0);
  for (std::size_t i = 0; i < n; ++i) {
    const auto cx = static_cast<std::size_t>((Y[2 * i] - xmin) / radius);
    const auto cy = static_cast<std::size_t>((Y[2 * i + 1] - ymin) / radius);
    cell[i] = cx + nx * cy;
    ++offsets[cell[i] + 1];
  }
  std::partial_sum(offsets.begin(), offsets.end(), offsets.begin());
  auto next = offsets;
  std::vector<std::size_t> order(n);
  for (std::size_t i = 0; i < n; ++i) {
    order[next[cell[i]]++] = i;
  }

  std::vector<double> sums(smallvis::worker_count(n, n_threads), 0.0);
  smallvis::parallel_for(
      n, n_threads, [&](std::size_t start, std::size_t end, std::size_t tid) {
        double total = 0.0;
        for (std::size_t i = start; i < end; ++i) {
          const auto cx = cell[i] % nx;
          const auto cy = cell[i] / nx;
          for (auto y = cy == 0 ? 0 : cy - 1; y < std::min(ny, cy + 2); ++y) {
            for (auto x = cx == 0 ? 0 : cx - 1; x < std::min(nx, cx + 2); ++x) {
              const auto bucket = x + nx * y;
              for (auto k = offsets[bucket]; k < offsets[bucket + 1]; ++k) {
                const auto j = order[k];
                if (i == j) {
                  continue;
                }
                const double dx = Y[2 * i] - Y[2 * j];
                const double dy = Y[2 * i + 1] - Y[2 * j + 1];
                const double d2 = dx * dx + dy * dy;
                if (d2 >= split.cutoff2) {
                  continue;
                }
                const double correction = split.correction(d2);
                if (split.gradient) {
                  gradient[2 * i] += correction * dx;
                  gradient[2 * i + 1] += correction * dy;
                } else {
                  total += correction;
                }
              }
            }
          }
        }
        sums[tid] = total;
      });
  return std::accumulate(sums.begin(), sums.end(), 0.0);
}

} // namespace

// [[Rcpp::export]]
SEXP bh_mixed_kl_cpp(const IntegerVector &indices, const IntegerVector &indptr,
                     const NumericVector &log_ratio,
                     const NumericVector &p_excess,
                     const NumericMatrix &embedding, double theta,
                     std::size_t n_threads, double eps, double lambda,
                     double p_mass, double plogp, bool gaussian,
                     bool gradient) {
  auto evaluate = [&](auto kernel) {
    if (gradient) {
      return bh_mixed_kl<true>(indices, indptr, log_ratio, p_excess, embedding,
                               theta, n_threads, eps, lambda, p_mass, plogp,
                               kernel);
    }
    return bh_mixed_kl<false>(indices, indptr, log_ratio, p_excess, embedding,
                              theta, n_threads, eps, lambda, p_mass, plogp,
                              kernel);
  };
  if (gaussian) {
    return evaluate(GaussianKernel{});
  }
  return evaluate(CauchyKernel{});
}

// FFT provides the kernel sums, including the log-weighted terms.
// [[Rcpp::export(rng = false)]]
SEXP fft_mixed_kl_cpp(const IntegerVector &indices, const IntegerVector &indptr,
                      const NumericVector &log_ratio,
                      const NumericVector &p_excess,
                      const NumericMatrix &embedding, double Z,
                      const NumericVector &log_mass,
                      Nullable<NumericMatrix> force,
                      Nullable<NumericMatrix> log_force, std::size_t n_threads,
                      double eps, double lambda, double p_mass, double plogp,
                      bool gaussian) {
  const auto Y = as<std::vector<double>>(transpose(embedding));
  const std::size_t n = Y.size() / 2;
  std::vector<KlRow> rows(n);
  double center[2] = {};
  for (std::size_t i = 0; i < n; ++i) {
    rows[i].log_mass = log_mass[i];
    center[0] += Y[2 * i];
    center[1] += Y[2 * i + 1];
  }
  center[0] /= n;
  center[1] /= n;
  auto evaluate = [&](auto kernel) {
    if (force.isNull()) {
      return mixed_kl_from_sums<false>(indices, indptr, log_ratio, p_excess, Y,
                                       center, rows, Z, n_threads, eps, lambda,
                                       p_mass, plogp, kernel);
    }
    const NumericMatrix raw_force(force);
    const NumericMatrix raw_log_force(log_force);
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t axis = 0; axis < 2; ++axis) {
        rows[i].force[axis] = -raw_force(i, axis);
        rows[i].log_force[axis] = -raw_log_force(i, axis);
      }
    }
    return mixed_kl_from_sums<true>(indices, indptr, log_ratio, p_excess, Y,
                                    center, rows, Z, n_threads, eps, lambda,
                                    p_mass, plogp, kernel);
  };
  if (gaussian) {
    return evaluate(GaussianKernel{});
  }
  return evaluate(CauchyKernel{});
}

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

// The FFT method supplies the slots of its validated sparse affinity matrix.
// [[Rcpp::export(rng = false)]]
NumericMatrix sne_sparse_attraction_cpp(const IntegerVector &indices,
                                        const IntegerVector &indptr,
                                        const NumericVector &P_data,
                                        const NumericMatrix &embedding,
                                        std::size_t n_threads,
                                        double alpha = 1.0) {
  if (embedding.ncol() != 2) {
    stop("FFT requires two output dimensions");
  }
  const auto Y = as<std::vector<double>>(transpose(embedding));
  std::vector<double> gradient(Y.size(), 0.0);
  auto accumulate = [&](const auto &kernel) {
    sparse_attraction(indices.begin(), indptr.begin(), P_data.begin(), Y,
                      n_threads, kernel, gradient);
  };
  if (alpha == 0.0) {
    accumulate(GaussianKernel{});
  } else if (alpha == 1.0) {
    accumulate(CauchyKernel{});
  } else {
    accumulate(HeavyTailKernel(alpha));
  }
  NumericMatrix result(2, embedding.nrow(), gradient.begin());
  return transpose(result);
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
  if (normalize && alpha == 0.0) {
    return bh_gradient<true>(indices, indptr, P_data, embedding, theta, 0.0,
                             n_threads, GaussianKernel{}, 0.0);
  }
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
  if (normalize && alpha == 0.0) {
    return bh_score(indices, indptr, P_data, embedding, theta, 0.0, n_threads,
                    GaussianKernel{}, true, 0.0);
  }
  if (normalize && alpha != 1.0) {
    return bh_score(indices, indptr, P_data, embedding, theta, eps, n_threads,
                    HeavyTailKernel(alpha), true, 0.0);
  }
  return bh_score(indices, indptr, P_data, embedding, theta, eps, n_threads,
                  CauchyKernel{}, normalize, lambda);
}

// [[Rcpp::export(rng = false)]]
NumericMatrix fft_bernoulli_grid_cpp(const NumericMatrix &d2,
                                     const std::string &method, double eps,
                                     double gr_eps, double a, double b,
                                     double radius, bool gradient,
                                     std::size_t n_threads) {
  NumericMatrix values(d2.nrow(), d2.ncol());
  const double *distance = d2.begin();
  double *output = values.begin();
  with_bernoulli_kernel(method, eps, gr_eps, a, b, [&](const auto &kernel) {
    const BernoulliGridKernel split(kernel, radius, gradient);
    smallvis::parallel_for(
        d2.size(), n_threads,
        [&](std::size_t start, std::size_t end, std::size_t) {
          for (std::size_t k = start; k < end; ++k) {
            output[k] = split.smooth(distance[k]);
          }
        });
    return 0;
  });
  return values;
}

// [[Rcpp::export(rng = false)]]
SEXP fft_bernoulli_finish_cpp(
    const IntegerVector &indices, const IntegerVector &indptr,
    const NumericVector &P_data, const NumericMatrix &embedding,
    Nullable<NumericMatrix> grid_force, double grid_mass,
    const std::string &method, double eps, double gr_eps, double a, double b,
    double gamma, double radius, std::size_t n_threads) {
  const auto Y = as<std::vector<double>>(transpose(embedding));
  return with_bernoulli_kernel(
      method, eps, gr_eps, a, b, [&](const auto &kernel) -> SEXP {
        const bool gradient = grid_force.isNotNull();
        const BernoulliGridKernel split(kernel, radius, gradient);
        if (gradient) {
          const NumericMatrix force(grid_force);
          auto result = as<std::vector<double>>(transpose(force));
          bernoulli_near_pairs(Y, radius, n_threads, split, result);
          for (auto &value : result) {
            value *= -gamma;
          }
          sparse_attraction(indices.begin(), indptr.begin(), P_data.begin(), Y,
                            n_threads, kernel, result);
          NumericMatrix gradient(2, embedding.nrow(), result.begin());
          return transpose(gradient);
        }

        std::vector<double> unused;
        const double near_mass =
            bernoulli_near_pairs(Y, radius, n_threads, split, unused);
        const double diag = method == "largevis"
                                ? embedding.nrow() * gamma * std::log1p(eps)
                                : 0.0;
        return wrap(sparse_log_weight_sum(indices.begin(), indptr.begin(),
                                          P_data.begin(), Y, n_threads, kernel,
                                          0.0,
                                          -gamma * (grid_mass + near_mass)) +
                    diag);
      });
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

// A missing force requests only the objective's sparse edge term.
// [[Rcpp::export(rng = false)]]
SEXP fft_divergence_cpp(const IntegerVector &indices,
                        const IntegerVector &indptr,
                        const NumericVector &edge_data,
                        const NumericMatrix &embedding, double Z,
                        Nullable<NumericMatrix> force, std::size_t n_threads,
                        const std::string &method) {
  const auto Y = as<std::vector<double>>(transpose(embedding));
  auto evaluate = [&](auto kind) -> SEXP {
    constexpr auto Kind = decltype(kind)::value;
    if (force.isNull()) {
      return wrap(divergence_score_from_mass<Kind>(indices, indptr, edge_data,
                                                   Y, Z, n_threads));
    }
    const NumericMatrix raw_force(force);
    const auto repulsion = as<std::vector<double>>(transpose(raw_force));
    // The FFT force uses query minus source, opposite to the tree sum.
    return divergence_gradient_from_sums<Kind>(indices, indptr, edge_data, Y, Z,
                                               repulsion, -1.0, n_threads);
  };
  if (method == "jssne") {
    return evaluate(
        std::integral_constant<Divergence, Divergence::JensenShannon>{});
  }
  if (method == "hlsne") {
    return evaluate(
        std::integral_constant<Divergence, Divergence::Hellinger>{});
  }
  if (method == "chsne") {
    return evaluate(
        std::integral_constant<Divergence, Divergence::ChiSquared>{});
  }
  stop("Unknown FFT divergence");
}

// [[Rcpp::export]]
NumericMatrix bh_jssne_gradient_cpp(const IntegerVector &indices,
                                    const IntegerVector &indptr,
                                    const NumericVector &P_data,
                                    const NumericMatrix &embedding,
                                    double theta, int n_threads) {
  return bh_divergence_gradient<Divergence::JensenShannon>(
      indices, indptr, P_data, embedding, theta, n_threads);
}

// [[Rcpp::export]]
double bh_jssne_score_cpp(const IntegerVector &indices,
                          const IntegerVector &indptr,
                          const NumericVector &P_data,
                          const NumericMatrix &embedding, double theta,
                          int n_threads) {
  return bh_divergence_score<Divergence::JensenShannon>(
      indices, indptr, P_data, embedding, theta, n_threads);
}

// [[Rcpp::export]]
NumericMatrix bh_hlsne_gradient_cpp(const IntegerVector &indices,
                                    const IntegerVector &indptr,
                                    const NumericVector &sqrt_P_data,
                                    const NumericMatrix &embedding,
                                    double theta, int n_threads) {
  return bh_divergence_gradient<Divergence::Hellinger>(
      indices, indptr, sqrt_P_data, embedding, theta, n_threads);
}

// [[Rcpp::export]]
double bh_hlsne_similarity_cpp(const IntegerVector &indices,
                               const IntegerVector &indptr,
                               const NumericVector &sqrt_P_data,
                               const NumericMatrix &embedding, double theta,
                               int n_threads) {
  return bh_divergence_score<Divergence::Hellinger>(
      indices, indptr, sqrt_P_data, embedding, theta, n_threads);
}

// [[Rcpp::export]]
NumericMatrix bh_chsne_gradient_cpp(const IntegerVector &indices,
                                    const IntegerVector &indptr,
                                    const NumericVector &squared_P_data,
                                    const NumericMatrix &embedding,
                                    double theta, int n_threads) {
  return bh_divergence_gradient<Divergence::ChiSquared>(
      indices, indptr, squared_P_data, embedding, theta, n_threads);
}

// [[Rcpp::export]]
double bh_chsne_score_cpp(const IntegerVector &indices,
                          const IntegerVector &indptr,
                          const NumericVector &squared_P_data,
                          const NumericMatrix &embedding, double theta,
                          int n_threads) {
  return bh_divergence_score<Divergence::ChiSquared>(
      indices, indptr, squared_P_data, embedding, theta, n_threads);
}
