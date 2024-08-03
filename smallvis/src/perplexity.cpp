#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <thread>
#include <vector>

#include <Rcpp.h>

using namespace Rcpp;

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
    beta[i] = (0.5 * perplexity * n) / sum_d2i;

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

    if (fabs(Hdiff) > tol) {
      bad_perp++;
      std::fill(W.begin() + idx, W.begin() + idx + n, 0.0);
      std::vector<int> knn_idx(n);
      std::iota(knn_idx.begin(), knn_idx.end(), 0);
      std::partial_sort(
          knn_idx.begin(),
          knn_idx.begin() + std::max(static_cast<int>(floor(perplexity)), 1),
          knn_idx.end(), [&D2i](int a, int b) { return D2i[a] < D2i[b]; });
      for (int k = 0; k < std::max(static_cast<int>(floor(perplexity)), 1);
           ++k) {
        W[idx + knn_idx[k]] = 1.0;
      }
    }
  }
}

// [[Rcpp::export]]
List find_beta_cpp(const NumericMatrix &X, double perplexity = 15,
                   double tol = 1e-5, int max_tries = 50,
                   std::size_t n_threads = 1) {

  std::size_t n = X.nrow();
  std::size_t d = X.ncol();
  double logU = log(perplexity);

  std::vector<double> X_vec(n * d);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < d; ++j) {
      X_vec[i * d + j] = X(i, j);
    }
  }

  std::vector<double> W(n * n, 0.0);
  std::vector<double> beta(n, 0.0);
  int bad_perp = 0;
  if (n_threads > 1) {
    std::size_t chunk_size = (n + n_threads - 1) / n_threads;
    std::vector<std::thread> threads;
    for (std::size_t t = 0; t < n_threads; ++t) {
      std::size_t start_row = t * chunk_size;
      std::size_t end_row = std::min(start_row + chunk_size, n);
      threads.emplace_back(find_beta, std::cref(X_vec), n, d, perplexity, logU,
                           tol, max_tries, std::ref(W), std::ref(beta),
                           std::ref(bad_perp), start_row, end_row);
    }
    for (auto &thread : threads) {
      thread.join();
    }
  } else {
    find_beta(X_vec, n, d, perplexity, logU, tol, max_tries, W, beta, bad_perp,
              0, n);
  }

  NumericMatrix W_mat(n, n);
  std::copy(W.begin(), W.end(), W_mat.begin());
  return List::create(Named("W") = transpose(W_mat), Named("beta") = beta,
                      Named("bad_perp") = bad_perp);
}