#include "threads.h"
#include <algorithm>
#include <cmath>
#include <limits>
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

// The caller excludes self by passing its index, or size for a kNN row.
double shifted_entropy(const double *distances, std::size_t size,
                       std::size_t self, double nearest, double beta) {
  double sum_w = 0.0;
  double sum_w_log_w = 0.0;
  for (std::size_t j = 0; j < size; ++j) {
    if (j != self) {
      const double log_w = -beta * (distances[j] - nearest);
      const double w = std::exp(log_w);
      sum_w += w;
      if (w > 0.0) {
        sum_w_log_w += w * log_w;
      }
    }
  }
  return std::log(sum_w) - sum_w_log_w / sum_w;
}

bool calibrate_row(const double *distances, std::size_t size, std::size_t self,
                   double nearest, double logU, double tol, int max_tries,
                   double &beta, double *row_scale, double *weights) {
  double lower = 0.0;
  double upper = std::numeric_limits<double>::infinity();
  double diff = shifted_entropy(distances, size, self, nearest, beta) - logU;

  for (int tries = 0;
       std::isfinite(diff) && std::abs(diff) > tol && tries < max_tries;
       ++tries) {
    if (diff > 0.0) {
      lower = beta;
      beta = std::isinf(upper) ? 2 * beta : (beta + upper) / 2;
    } else {
      upper = beta;
      beta = (beta + lower) / 2;
    }
    diff = shifted_entropy(distances, size, self, nearest, beta) - logU;
  }

  if (!std::isfinite(diff) || std::abs(diff) > tol) {
    return false;
  }
  const double raw_scale = std::exp(-beta * nearest);
  if (row_scale == nullptr && raw_scale == 0.0) {
    return false;
  }
  if (row_scale != nullptr) {
    *row_scale = raw_scale;
  }
  const double offset = row_scale == nullptr ? 0.0 : nearest;
  for (std::size_t j = 0; j < size; ++j) {
    weights[j] = j == self ? 0.0 : std::exp(-beta * (distances[j] - offset));
  }
  return true;
}

std::vector<int>
fallback_rows(const std::vector<unsigned char> &used_fallback) {
  std::vector<int> rows;
  for (std::size_t i = 0; i < used_fallback.size(); ++i) {
    if (used_fallback[i]) {
      rows.push_back(static_cast<int>(i + 1));
    }
  }
  return rows;
}

void find_beta_knn(std::vector<double> &knn_distances, std::size_t k,
                   double perplexity, double logU, double tol, int max_tries,
                   std::vector<double> &W, std::vector<double> &beta,
                   double *row_scale, std::vector<unsigned char> &used_fallback,
                   std::size_t start_row, std::size_t end_row) {
  for (std::size_t i = start_row; i < end_row; ++i) {
    const std::size_t idx = i * k;
    auto *distances = knn_distances.data() + idx;
    double sum_d2i = 0.0;
    double nearest = std::numeric_limits<double>::infinity();
    for (std::size_t j = 0; j < k; ++j) {
      distances[j] *= distances[j];
      sum_d2i += distances[j];
      nearest = std::min(nearest, distances[j]);
    }

    beta[i] = sum_d2i == 0.0 ? 0.0 : (0.5 * perplexity * k) / sum_d2i;
    if (!calibrate_row(distances, k, k, nearest, logU, tol, max_tries, beta[i],
                       row_scale == nullptr ? nullptr : row_scale + i,
                       W.data() + idx)) {
      used_fallback[i] = 1;
      std::fill_n(W.data() + idx, k, 1.0 / k);
    }
  }
}

void find_beta(const std::vector<double> &data, std::size_t n, std::size_t d,
               double perplexity, double logU, double tol, int max_tries,
               double *W, std::vector<double> &beta, double *row_scale,
               std::vector<unsigned char> &used_fallback, std::size_t start_row,
               std::size_t end_row) {

  // Trials use contiguous per-worker storage. Scatter each completed row into
  // R's column-major output once, avoiding quadratic staging and transposes.
  std::vector<double> D2i(n);
  std::vector<double> Wi(n);
  for (std::size_t i = start_row; i < end_row; ++i) {
    const std::size_t i_d = i * d;

    double sum_d2i = 0.0;
    double nearest = std::numeric_limits<double>::infinity();

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
      if (i != j) {
        nearest = std::min(nearest, dist);
      }
    }

    // Initial guess for beta: 0.5 * perplexity / mean(D2i)
    beta[i] = sum_d2i == 0.0 ? 0.0 : (0.5 * perplexity * n) / sum_d2i;
    if (!calibrate_row(D2i.data(), n, i, nearest, logU, tol, max_tries, beta[i],
                       row_scale == nullptr ? nullptr : row_scale + i,
                       Wi.data())) {
      used_fallback[i] = 1;
      std::fill(Wi.begin(), Wi.end(), 0.0);
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
        Wi[knn_idx[k]] = 1.0;
      }
    }
    for (std::size_t j = 0; j < n; ++j) {
      W[i + j * n] = Wi[j];
    }
  }
}

// [[Rcpp::export]]
List find_beta_knn_cpp(const NumericMatrix &knn_distances,
                       const IntegerMatrix &knn_indices, double perplexity = 15,
                       double tol = 1e-5, int max_tries = 50,
                       bool ret_sparse = false, std::size_t n_threads = 1,
                       bool shift_rows = false) {

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
  std::vector<double> row_scale(shift_rows ? n : 0, 1.0);
  double *scale = shift_rows ? row_scale.data() : nullptr;
  // Each worker owns a disjoint range of row flags.
  std::vector<unsigned char> used_fallback(n, 0);
  auto worker = [&](std::size_t start, std::size_t end, std::size_t) {
    find_beta_knn(knn_distances_vec, k, perplexity, logU, tol, max_tries, W,
                  beta, scale, used_fallback, start, end);
  };
  smallvis::parallel_for(n, n_threads, worker);
  const auto rows = fallback_rows(used_fallback);
  const auto bad_perp = static_cast<int>(rows.size());

  if (ret_sparse) {
    NumericVector P(n * k);
    std::copy(W.begin(), W.end(), P.begin());
    return List::create(
        Named("P") = P, Named("beta") = beta, Named("bad_perp") = bad_perp,
        Named("fallback_rows") = rows, Named("row_scale") = row_scale);
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
    return List::create(
        Named("P") = P, Named("beta") = beta, Named("bad_perp") = bad_perp,
        Named("fallback_rows") = rows, Named("row_scale") = row_scale);
  }
}

// [[Rcpp::export]]
List find_beta_cpp(const NumericMatrix &X, double perplexity = 15,
                   double tol = 1e-5, int max_tries = 50,
                   std::size_t n_threads = 1, bool shift_rows = false) {

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

  NumericMatrix W(n, n);
  double *W_data = W.begin();
  std::vector<double> beta(n, 0.0);
  std::vector<double> row_scale(shift_rows ? n : 0, 1.0);
  double *scale = shift_rows ? row_scale.data() : nullptr;
  std::vector<unsigned char> used_fallback(n, 0);
  auto worker = [&](std::size_t start, std::size_t end, std::size_t) {
    find_beta(X_vec, n, d, perplexity, logU, tol, max_tries, W_data, beta,
              scale, used_fallback, start, end);
  };
  smallvis::parallel_for(n, n_threads, worker);
  const auto rows = fallback_rows(used_fallback);
  const auto bad_perp = static_cast<int>(rows.size());

  return List::create(
      Named("W") = W, Named("beta") = beta, Named("bad_perp") = bad_perp,
      Named("fallback_rows") = rows, Named("row_scale") = row_scale);
}
