#include "bh.h"
#include <Rcpp.h>

using namespace Rcpp;

void _bh_tsne_negative_gradient_single(const smallvis::QuadTree::Node &node,
                                       double point_x, double point_y,
                                       double theta2, double eps,
                                       double &gradient_x, double &gradient_y,
                                       double &Zi) {
  // Ensure we do not process empty nodes or simple self-interactions
  if (node.is_skippable(point_x, point_y, eps)) {
    return;
  }

  // Compute the squared Euclidean distance between the point and the center of
  // mass
  double diff_x = node.center_of_mass_x - point_x;
  double diff_y = node.center_of_mass_y - point_y;
  double d2_ij = eps + diff_x * diff_x + diff_y * diff_y;

  // Check if we can use this node as a summary
  if (node.is_summary(d2_ij, theta2)) {
    double w_ij = 1.0 / (1.0 + d2_ij);
    double nw_ij = node.num_points * w_ij;
    Zi += nw_ij;

    nw_ij *= w_ij;
    // un-normalized -ve grad W*W*dY
    gradient_x += nw_ij * diff_x;
    gradient_y += nw_ij * diff_y;

    return;
  }

  // Recursively apply Barnes-Hut to the children
  for (const auto &child : node.children) {
    _bh_tsne_negative_gradient_single(*child, point_x, point_y, theta2, eps,
                                      gradient_x, gradient_y, Zi);
  }
}

// Function to estimate the negative gradient using the Barnes-Hut approximation
void bh_tsne_negative_gradient(const smallvis::QuadTree &tree,
                               const std::vector<double> &embedding,
                               double theta2, double eps, std::size_t n_threads,
                               std::vector<double> &gradient) {
  std::size_t embedding_len = embedding.size();
  std::vector<double> Zi(std::max(static_cast<std::size_t>(1), n_threads), 0.0);

  const smallvis::QuadTree::Node *root = tree.root.get();

  // Function to calculate the negative gradient for a single point
  auto worker = [&](std::size_t start, std::size_t end, std::size_t thread_id) {
    for (std::size_t i = start * 2; i < end * 2; i += 2) {
      _bh_tsne_negative_gradient_single(*root, embedding[i], embedding[i + 1],
                                        theta2, eps, gradient[i],
                                        gradient[i + 1], Zi[thread_id]);
    }
  };
  smallvis::parallel_for(embedding_len / 2, n_threads, worker);

  double Z = std::accumulate(Zi.begin(), Zi.end(), eps);
  // Normalize the gradient
  for (std::size_t i = 0; i < embedding_len; i += 2) {
    gradient[i] /= Z;
    gradient[i + 1] /= Z;
  }
}

void bh_tsne_positive_gradient(const std::vector<std::size_t> &indices,
                               const std::vector<std::size_t> &indptr,
                               const std::vector<double> &P_data,
                               const std::vector<double> &embedding,
                               std::size_t n_threads,
                               std::vector<double> &gradient) {

  auto worker = [&](std::size_t start, std::size_t end) {
    for (std::size_t i = start, i2 = start * 2; i < end; ++i, i2 += 2) {
      for (std::size_t k = indptr[i]; k < indptr[i + 1]; ++k) {
        std::size_t j2 = indices[k] * 2;

        // Compute the direction of the points' attraction and the squared
        // Euclidean distance
        double diff_x = embedding[i2] - embedding[j2];
        double diff_y = embedding[i2 + 1] - embedding[j2 + 1];
        double w_ij_p_ij =
            P_data[k] / (1.0 + diff_x * diff_x + diff_y * diff_y);

        // Compute F_{attr} of point `j` on point `i`
        // W x P x dY
        gradient[i2] += w_ij_p_ij * diff_x;
        gradient[i2 + 1] += w_ij_p_ij * diff_y;
      }
    }
  };
  smallvis::parallel_for(gradient.size() / 2, n_threads, worker);
}

double _bh_Zi(const smallvis::QuadTree::Node &node, double point_x,
              double point_y, double theta2, double eps) {
  // Ensure we do not process empty nodes or simple self-interactions
  if (node.is_skippable(point_x, point_y, eps)) {
    return 0.0;
  }

  double diff_x = node.center_of_mass_x - point_x;
  double diff_y = node.center_of_mass_y - point_y;
  double d2_ij = eps + diff_x * diff_x + diff_y * diff_y;

  if (node.is_summary(d2_ij, theta2)) {
    return node.num_points / (1.0 + d2_ij);
  }

  // Recursively apply Barnes-Hut to the children
  double Zi = 0.0;
  for (const auto &child : node.children) {
    Zi += _bh_Zi(*child, point_x, point_y, theta2, eps);
  }
  return Zi;
}

double bh_Z(const smallvis::QuadTree &tree,
            const std::vector<double> &embedding, double theta2, double eps,
            std::size_t n_threads) {

  std::vector<double> Zi(std::max(static_cast<std::size_t>(1), n_threads), 0.0);
  const smallvis::QuadTree::Node *root = tree.root.get();
  auto worker = [&](std::size_t start, std::size_t end, std::size_t thread_id) {
    const std::size_t end2 = end * 2;
    for (std::size_t i = start * 2; i < end2; i += 2) {
      Zi[thread_id] +=
          _bh_Zi(*root, embedding[i], embedding[i + 1], theta2, eps);
    }
  };
  smallvis::parallel_for(embedding.size() / 2, n_threads, worker);

  return std::accumulate(Zi.begin(), Zi.end(), eps);
}

double bh_plogq(const std::vector<std::size_t> &indices,
                const std::vector<std::size_t> &indptr,
                const std::vector<double> &P_data,
                const std::vector<double> &embedding, double theta, double eps,
                std::size_t n_threads) {
  smallvis::QuadTree tree(embedding, eps);
  double mlogZ =
      -std::log(bh_Z(tree, embedding, theta * theta, eps, n_threads));

  std::vector<double> partial_plogq(
      std::max(static_cast<std::size_t>(1), n_threads));
  auto worker = [&](std::size_t start, std::size_t end, std::size_t thread_id) {
    double window_plogq = 0.0;
    for (std::size_t i = start, i2 = start * 2; i < end; ++i, i2 += 2) {
      for (std::size_t k = indptr[i]; k < indptr[i + 1]; ++k) {
        std::size_t j2 = indices[k] * 2;

        // p * log(q) = p * log(w/Z) = p * (log w - log Z) =
        // p * (log(1/(1+d2)) - log Z) = p * (-log(1+d2) - log Z)
        double diff_x = embedding[i2] - embedding[j2];
        double diff_y = embedding[i2 + 1] - embedding[j2 + 1];

        window_plogq +=
            P_data[k] *
            (mlogZ - std::log(1.0 + diff_x * diff_x + diff_y * diff_y));
      }
    }
    partial_plogq[thread_id] += window_plogq;
  };
  smallvis::parallel_for(embedding.size() / 2, n_threads, worker);

  return std::accumulate(partial_plogq.begin(), partial_plogq.end(), 0.0);
}

std::vector<double> bh_tsne_gradient(const std::vector<std::size_t> &indices,
                                     const std::vector<std::size_t> &indptr,
                                     const std::vector<double> &P_data,
                                     const std::vector<double> &embedding,
                                     double theta, double eps,
                                     std::size_t n_threads) {

  smallvis::QuadTree tree(embedding, eps);
  std::vector<double> gradient(embedding.size(), 0.0);

  // don't re-arrange the order of the following two function calls!
  bh_tsne_negative_gradient(tree, embedding, theta * theta, eps, n_threads,
                            gradient);
  bh_tsne_positive_gradient(indices, indptr, P_data, embedding, n_threads,
                            gradient);
  return gradient;
}

// [[Rcpp::export]]
NumericMatrix bh_tsne_gradient_cpp(const IntegerVector &indices,
                                   const IntegerVector &indptr,
                                   const NumericVector &P_data,
                                   const NumericMatrix &embedding,
                                   double theta = 0.5, double eps = 1e-16,
                                   std::size_t n_threads = 1) {
  auto indices_cpp = as<std::vector<std::size_t>>(indices);
  auto indptr_cpp = as<std::vector<std::size_t>>(indptr);
  auto P_data_cpp = as<std::vector<double>>(P_data);

  auto embedding_vec = as<std::vector<double>>(transpose(embedding));

  std::vector<double> gradient_cpp =
      bh_tsne_gradient(indices_cpp, indptr_cpp, P_data_cpp, embedding_vec,
                       theta, eps, n_threads);

  NumericMatrix gradient(2, embedding.nrow(), gradient_cpp.begin());
  return transpose(gradient);
}

// [[Rcpp::export]]
double bh_plogq_cpp(const IntegerVector &indices, const IntegerVector &indptr,
                    const NumericVector &P_data, const NumericMatrix &embedding,
                    double theta = 0.5, double eps = 1e-16,
                    std::size_t n_threads = 1) {
  auto indices_cpp = as<std::vector<std::size_t>>(indices);
  auto indptr_cpp = as<std::vector<std::size_t>>(indptr);
  auto P_data_cpp = as<std::vector<double>>(P_data);

  auto embedding_vec = as<std::vector<double>>(transpose(embedding));

  return bh_plogq(indices_cpp, indptr_cpp, P_data_cpp, embedding_vec, theta,
                  eps, n_threads);
}
