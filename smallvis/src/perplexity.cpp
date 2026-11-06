#include "threads.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

#include <Rcpp.h>

using namespace Rcpp;

void validate_perplexity(double perplexity, std::size_t available) {
  if (!std::isfinite(perplexity) || perplexity < 1.0 ||
      perplexity > static_cast<double>(available)) {
    stop("Perplexity must be finite and between 1 and the number of "
         "non-self neighbors");
  }
}

void find_beta_knn(std::vector<double> &knn_distances, std::size_t k,
                   double perplexity, double logU, double tol, int max_tries,
                   std::vector<double> &W, std::vector<double> &beta,
                   int &bad_perp, std::size_t start_row, std::size_t end_row) {

  for (std::size_t i = start_row; i < end_row; ++i) {
    const std::size_t idx = i * k;
    double betamin = -std::numeric_limits<double>::infinity();
    double betamax = std::numeric_limits<double>::infinity();

    // Square the distances
    for (std::size_t j = 0; j < k; ++j) {
      knn_distances[idx + j] *= knn_distances[idx + j];
    }

    // Initial guess for beta: 0.5 * perplexity / mean(knn_distances)
    double sum_d2i = std::accumulate(knn_distances.begin() + idx,
                                     knn_distances.begin() + idx + k, 0.0);
    beta[i] = sum_d2i == 0.0 ? 0.0 : (0.5 * perplexity * k) / sum_d2i;

    double Z = 0.0;
    double entropy = 0.0;
    for (std::size_t j = 0; j < k; ++j) {
      W[idx + j] = exp(-knn_distances[idx + j] * beta[i]);
      entropy += knn_distances[idx + j] * W[idx + j];
      Z += W[idx + j];
    }
    if (Z == 0.0) {
      entropy = 0.0;
    } else {
      entropy = (entropy / Z) * beta[i] + log(Z);
    }

    double Hdiff = entropy - logU;
    int tries = 0;
    while (fabs(Hdiff) > tol && tries < max_tries) {
      if (Hdiff > 0) {
        betamin = beta[i];
        if (std::isinf(betamax)) {
          beta[i] *= 2;
        } else {
          beta[i] = (beta[i] + betamax) / 2;
        }
      } else {
        betamax = beta[i];
        if (std::isinf(betamin)) {
          beta[i] /= 2;
        } else {
          beta[i] = (beta[i] + betamin) / 2;
        }
      }

      Z = 0.0;
      entropy = 0.0;
      for (std::size_t j = 0; j < k; ++j) {
        W[idx + j] = exp(-knn_distances[idx + j] * beta[i]);
        entropy += knn_distances[idx + j] * W[idx + j];
        Z += W[idx + j];
      }
      if (Z == 0.0) {
        entropy = 0.0;
      } else {
        entropy = (entropy / Z) * beta[i] + log(Z);
      }

      Hdiff = entropy - logU;
      tries++;
    }

    if (!std::isfinite(Hdiff) || fabs(Hdiff) > tol) {
      bad_perp++;
      std::fill(W.begin() + idx, W.begin() + idx + k, 1.0 / k);
    }
  }
}

void find_beta(const std::vector<double> &data, std::size_t n, std::size_t d,
               double perplexity, double logU, double tol, int max_tries,
               std::vector<double> &W, std::vector<double> &beta, int &bad_perp,
               std::size_t start_row, std::size_t end_row) {

  for (std::size_t i = start_row; i < end_row; ++i) {
    const std::size_t i_d = i * d;
    const std::size_t idx = i * n;

    std::vector<double> D2i(n);
    double sum_d2i = 0.0;

    // D2
    for (std::size_t j = 0; j < n; ++j) {
      const std::size_t j_d = j * d;
      double dist = 0.0;
      for (std::size_t k = 0; k < d; ++k) {
        double diff = data[i_d + k] - data[j_d + k];
        dist += diff * diff;
      }
      sum_d2i += dist;
      D2i[j] = dist;
    }

    double betamin = -std::numeric_limits<double>::infinity();
    double betamax = std::numeric_limits<double>::infinity();

    // Initial guess for beta: 0.5 * perplexity / mean(D2i)
    beta[i] = sum_d2i == 0.0 ? 0.0 : (0.5 * perplexity * n) / sum_d2i;

    double Z = 0.0;
    double entropy = 0.0;
    for (std::size_t j = 0; j < n; ++j) {
      if (i == j) {
        W[idx + j] = 0.0;
        continue;
      }
      W[idx + j] = exp(-D2i[j] * beta[i]);
      entropy += D2i[j] * W[idx + j];
      Z += W[idx + j];
    }
    if (Z == 0.0) {
      entropy = 0.0;
    } else {
      entropy = (entropy / Z) * beta[i] + log(Z);
    }

    double Hdiff = entropy - logU;

    int tries = 0;
    while (fabs(Hdiff) > tol && tries < max_tries) {
      if (Hdiff > 0) {
        betamin = beta[i];
        if (std::isinf(betamax)) {
          beta[i] *= 2;
        } else {
          beta[i] = (beta[i] + betamax) / 2;
        }
      } else {
        betamax = beta[i];
        if (std::isinf(betamin)) {
          beta[i] /= 2;
        } else {
          beta[i] = (beta[i] + betamin) / 2;
        }
      }

      Z = 0.0;
      entropy = 0.0;
      for (std::size_t j = 0; j < n; ++j) {
        if (i == j) {
          W[idx + j] = 0.0;
          continue;
        }
        W[idx + j] = exp(-D2i[j] * beta[i]);
        entropy += D2i[j] * W[idx + j];
        Z += W[idx + j];
      }
      if (Z == 0.0) {
        entropy = 0.0;
      } else {
        entropy = (entropy / Z) * beta[i] + log(Z);
      }

      Hdiff = entropy - logU;
      tries++;
    }

    if (!std::isfinite(Hdiff) || fabs(Hdiff) > tol) {
      bad_perp++;
      std::fill(W.begin() + idx, W.begin() + idx + n, 0.0);
      std::vector<std::size_t> knn_idx;
      knn_idx.reserve(n - 1);
      for (std::size_t j = 0; j < n; ++j) {
        if (j != i) {
          knn_idx.push_back(j);
        }
      }
      const auto n_neighbors =
          std::min(n - 1, static_cast<std::size_t>(std::floor(perplexity)));
      std::partial_sort(knn_idx.begin(), knn_idx.begin() + n_neighbors,
                        knn_idx.end(), [&D2i](std::size_t a, std::size_t b) {
                          return D2i[a] == D2i[b] ? a < b : D2i[a] < D2i[b];
                        });
      for (std::size_t k = 0; k < n_neighbors; ++k) {
        W[idx + knn_idx[k]] = 1.0;
      }
    }
  }
}

// [[Rcpp::export]]
List find_beta_knn_cpp(const NumericMatrix &knn_distances,
                       const IntegerMatrix &knn_indices, double perplexity = 15,
                       double tol = 1e-5, int max_tries = 50,
                       bool ret_sparse = false, std::size_t n_threads = 1) {

  const std::size_t n = knn_distances.nrow();
  const std::size_t k = knn_distances.ncol();
  validate_perplexity(perplexity, std::min(k, n > 0 ? n - 1 : 0));
  if (knn_indices.nrow() != knn_distances.nrow() ||
      knn_indices.ncol() != knn_distances.ncol()) {
    stop("Neighbor index and distance matrices must have the same dimensions");
  }
  const double logU = log(perplexity);

  // Flatten the knn_distances matrix into a vector
  std::vector<double> knn_distances_vec(n * k);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < k; ++j) {
      knn_distances_vec[i * k + j] = knn_distances(i, j);
    }
  }

  std::vector<double> W(n * k, 0.0);
  std::vector<double> beta(n, 0.0);
  std::vector<int> failures(smallvis::worker_count(n, n_threads), 0);
  auto worker = [&](std::size_t start, std::size_t end, std::size_t id) {
    find_beta_knn(knn_distances_vec, k, perplexity, logU, tol, max_tries, W,
                  beta, failures[id], start, end);
  };
  smallvis::parallel_for(n, n_threads, worker);
  const int bad_perp = std::accumulate(failures.begin(), failures.end(), 0);

  if (ret_sparse) {
    NumericVector P(n * k);
    std::copy(W.begin(), W.end(), P.begin());
    return List::create(Named("P") = P, Named("beta") = beta,
                        Named("bad_perp") = bad_perp);
  } else {
    NumericMatrix P(n, n);
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = 0; j < k; ++j) {
        const auto index = knn_indices(i, j);
        if (index < 1 || static_cast<std::size_t>(index) > n) {
          stop("Neighbor index is out of range");
        }
        P(i, index - 1) = W[i * k + j];
      }
    }
    return List::create(Named("P") = P, Named("beta") = beta,
                        Named("bad_perp") = bad_perp);
  }
}

// [[Rcpp::export]]
List find_beta_cpp(const NumericMatrix &X, double perplexity = 15,
                   double tol = 1e-5, int max_tries = 50,
                   std::size_t n_threads = 1) {

  std::size_t n = X.nrow();
  std::size_t d = X.ncol();
  validate_perplexity(perplexity, n > 0 ? n - 1 : 0);
  double logU = log(perplexity);

  std::vector<double> X_vec(n * d);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < d; ++j) {
      X_vec[i * d + j] = X(i, j);
    }
  }

  std::vector<double> W(n * n, 0.0);
  std::vector<double> beta(n, 0.0);
  std::vector<int> failures(smallvis::worker_count(n, n_threads), 0);
  auto worker = [&](std::size_t start, std::size_t end, std::size_t id) {
    find_beta(X_vec, n, d, perplexity, logU, tol, max_tries, W, beta,
              failures[id], start, end);
  };
  smallvis::parallel_for(n, n_threads, worker);
  const int bad_perp = std::accumulate(failures.begin(), failures.end(), 0);

  NumericMatrix W_mat(n, n);
  std::copy(W.begin(), W.end(), W_mat.begin());
  return List::create(Named("W") = transpose(W_mat), Named("beta") = beta,
                      Named("bad_perp") = bad_perp);
}
