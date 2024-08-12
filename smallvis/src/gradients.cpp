#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <thread>
#include <vector>

using namespace Rcpp;

void mmds_grad(const std::vector<double> &R, const std::vector<double> &D,
               const std::vector<double> &Y, double eps,
               std::vector<double> &gradient, std::size_t start,
               std::size_t end, std::size_t n, std::size_t d) {

  for (std::size_t i = start; i < end; ++i) {
    const std::size_t i_d = i * d;
    const std::size_t i_n = i * n;
    for (std::size_t j = 0; j < n; ++j) {
      if (i == j) {
        continue;
      }
      const std::size_t ij = i_n + j;
      double k_ij = 4.0 * (D[ij] - R[ij]) / (D[ij] + eps);
      for (std::size_t k = 0; k < d; ++k) {
        gradient[i_d + k] += k_ij * (Y[i_d + k] - Y[j * d + k]);
      }
    }
  }
}

void tsne_grad(const std::vector<double> &P, const std::vector<double> &W,
               double Z, const std::vector<double> &Y,
               std::vector<double> &gradient, std::size_t start,
               std::size_t end, std::size_t n) {

  const double Z4 = 4.0 / Z;

  for (std::size_t i = start; i < end; ++i) {
    const std::size_t i_d = i + i;
    const std::size_t i_n = i * n;
    for (std::size_t j = 0; j < n; ++j) {
      if (i == j) {
        continue;
      }
      const std::size_t ij = i_n + j;
      double k_ij = Z4 * W[ij] * (Z * P[ij] - W[ij]);
      gradient[i_d] += k_ij * (Y[i_d] - Y[j * 2]);
      gradient[i_d + 1] += k_ij * (Y[i_d + 1] - Y[j * 2 + 1]);
    }
  }
}

void tweight(const std::vector<double> &data, std::vector<double> &dist_matrix,
             std::size_t start_row, std::size_t end_row, std::size_t n,
             std::size_t d) {
  for (std::size_t i = start_row; i < end_row; ++i) {
    const std::size_t i_d = i * d;
    const std::size_t i_n = i * n;
    for (std::size_t j = 0; j < n; ++j) {
      const std::size_t j_d = j * d;
      double dist = 0.0;
      for (std::size_t k = 0; k < d; ++k) {
        double diff = data[i_d + k] - data[j_d + k];
        dist += diff * diff;
      }
      dist_matrix[i_n + j] = 1.0 / (1.0 + dist);
    }
  }
}

void d2(const std::vector<double> &data, std::vector<double> &dist_matrix,
        std::size_t start_row, std::size_t end_row, std::size_t n,
        std::size_t d) {
  for (std::size_t i = start_row; i < end_row; ++i) {
    const std::size_t i_d = i * d;
    const std::size_t i_n = i * n;
    for (std::size_t j = 0; j < n; ++j) {
      const std::size_t j_d = j * d;
      double dist = 0.0;
      for (std::size_t k = 0; k < d; ++k) {
        double diff = data[i_d + k] - data[j_d + k];
        dist += diff * diff;
      }
      dist_matrix[i_n + j] = dist;
    }
  }
}

void dist(const std::vector<double> &data, std::vector<double> &dist_matrix,
          std::size_t start_row, std::size_t end_row, std::size_t n,
          std::size_t d) {
  for (std::size_t i = start_row; i < end_row; ++i) {
    const std::size_t i_d = i * d;
    const std::size_t i_n = i * n;
    for (std::size_t j = 0; j < n; ++j) {
      const std::size_t j_d = j * d;
      double dist = 0.0;
      for (std::size_t k = 0; k < d; ++k) {
        double diff = data[i_d + k] - data[j_d + k];
        dist += diff * diff;
      }
      dist_matrix[i_n + j] = sqrt(dist);
    }
  }
}

void d2_to_tweight(const std::vector<double> &dist_matrix,
                   std::vector<double> &transformed_matrix, int start,
                   int end) {
  for (int idx = start; idx < end; ++idx) {
    transformed_matrix[idx] = 1.0 / (1.0 + dist_matrix[idx]);
  }
}

// [[Rcpp::export]]
NumericMatrix dist2_cpp(NumericMatrix input, std::size_t n_threads = 1) {
  if (n_threads == 0) {
    n_threads = 1;
  }
  std::size_t n = input.nrow();
  std::size_t d = input.ncol();

  std::vector<double> transposed_data(n * d);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < d; ++j) {
      transposed_data[i * d + j] = input(i, j);
    }
  }

  std::vector<double> dist_matrix(n * n, 0.0);

  std::size_t chunk_size = (n + n_threads - 1) / n_threads;

  std::vector<std::thread> threads;
  for (std::size_t t = 0; t < n_threads; ++t) {
    std::size_t start_row = t * chunk_size;
    std::size_t end_row = std::min(start_row + chunk_size, n);
    threads.emplace_back(d2, std::cref(transposed_data), std::ref(dist_matrix),
                         start_row, end_row, n, d);
  }
  for (auto &thread : threads) {
    thread.join();
  }

  NumericMatrix result(n, n);
  std::copy(dist_matrix.begin(), dist_matrix.end(), result.begin());

  return result;
}

// [[Rcpp::export]]
NumericMatrix dist_cpp(NumericMatrix input, std::size_t n_threads = 1) {
  if (n_threads == 0) {
    n_threads = 1;
  }
  std::size_t n = input.nrow();
  std::size_t d = input.ncol();

  std::vector<double> transposed_data(n * d);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < d; ++j) {
      transposed_data[i * d + j] = input(i, j);
    }
  }

  std::vector<double> dist_matrix(n * n, 0.0);

  std::size_t chunk_size = (n + n_threads - 1) / n_threads;

  std::vector<std::thread> threads;
  for (std::size_t t = 0; t < n_threads; ++t) {
    std::size_t start_row = t * chunk_size;
    std::size_t end_row = std::min(start_row + chunk_size, n);
    threads.emplace_back(dist, std::cref(transposed_data),
                         std::ref(dist_matrix), start_row, end_row, n, d);
  }

  for (auto &thread : threads) {
    thread.join();
  }

  NumericMatrix result(n, n);
  std::copy(dist_matrix.begin(), dist_matrix.end(), result.begin());

  return result;
}

// [[Rcpp::export]]
NumericMatrix tweight_cpp(NumericMatrix input, std::size_t n_threads = 1) {
  if (n_threads == 0) {
    n_threads = 1;
  }
  std::size_t n = input.nrow();
  std::size_t d = input.ncol();

  std::vector<double> transposed_data(n * d);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < d; ++j) {
      transposed_data[i * d + j] = input(i, j);
    }
  }

  std::vector<double> dist_matrix(n * n, 0.0);

  std::size_t chunk_size = (n + n_threads - 1) / n_threads;

  std::vector<std::thread> threads;
  for (std::size_t t = 0; t < n_threads; ++t) {
    std::size_t start_row = t * chunk_size;
    std::size_t end_row = std::min(start_row + chunk_size, n);
    threads.emplace_back(tweight, std::cref(transposed_data),
                         std::ref(dist_matrix), start_row, end_row, n, d);
  }

  for (auto &thread : threads) {
    thread.join();
  }

  NumericMatrix result(n, n);
  std::copy(dist_matrix.begin(), dist_matrix.end(), result.begin());

  return result;
}

// [[Rcpp::export]]
NumericMatrix d2_to_tweight_cpp(NumericMatrix dist_matrix, int n_threads) {
  if (n_threads == 0) {
    n_threads = 1;
  }
  int n = dist_matrix.nrow();
  int total_elements = n * n;

  std::vector<double> dist_matrix_vec(dist_matrix.begin(), dist_matrix.end());

  std::vector<double> transformed_matrix(total_elements, 0.0);

  int chunk_size = (total_elements + n_threads - 1) / n_threads;

  std::vector<std::thread> threads;
  for (int t = 0; t < n_threads; ++t) {
    int start = t * chunk_size;
    int end = std::min(start + chunk_size, total_elements);
    threads.emplace_back(d2_to_tweight, std::cref(dist_matrix_vec),
                         std::ref(transformed_matrix), start, end);
  }

  for (auto &thread : threads) {
    thread.join();
  }

  NumericMatrix result(n, n);
  std::copy(transformed_matrix.begin(), transformed_matrix.end(),
            result.begin());

  return result;
}

// [[Rcpp::export]]
NumericMatrix tsne_grad_cpp(const NumericMatrix &P, const NumericMatrix &W,
                            double Z, const NumericMatrix &Y,
                            std::size_t n_threads) {
  std::size_t n = Y.nrow();

  std::vector<double> P_vec(P.begin(), P.end());
  std::vector<double> W_vec(W.begin(), W.end());

  std::vector<double> Y_vec(n * 2);
  for (std::size_t i = 0; i < n; ++i) {
    Y_vec[i * 2] = Y(i, 0);
    Y_vec[i * 2 + 1] = Y(i, 1);
  }

  std::vector<double> gradient_vec(n * 2, 0.0);

  if (n_threads > 1) {
    std::size_t chunk_size = (n + n_threads - 1) / n_threads;
    std::vector<std::thread> threads;
    for (std::size_t t = 0; t < n_threads; ++t) {
      std::size_t start_row = t * chunk_size;
      std::size_t end_row = std::min(start_row + chunk_size, n);
      threads.emplace_back(tsne_grad, std::cref(P_vec), std::cref(W_vec), Z,
                           std::cref(Y_vec), std::ref(gradient_vec), start_row,
                           end_row, n);
    }
    for (auto &thread : threads) {
      thread.join();
    }
  } else {
    tsne_grad(P_vec, W_vec, Z, Y_vec, gradient_vec, 0, n, n);
  }

  NumericMatrix gradient(n, 2);
  for (std::size_t i = 0; i < n; ++i) {
    gradient(i, 0) = gradient_vec[i * 2];
    gradient(i, 1) = gradient_vec[i * 2 + 1];
  }

  return gradient;
}

// [[Rcpp::export]]
NumericMatrix mmds_grad_cpp(const NumericMatrix &R, const NumericMatrix &D,
                            const NumericMatrix &Y, double eps,
                            std::size_t n_threads) {
  std::size_t n = Y.nrow();
  std::size_t d = Y.ncol();

  std::vector<double> R_vec(R.begin(), R.end());
  std::vector<double> D_vec(D.begin(), D.end());

  std::vector<double> Y_vec(n * d);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < d; ++j) {
      Y_vec[i * d + j] = Y(i, j);
    }
  }

  std::vector<double> gradient_vec(n * d, 0.0);

  if (n_threads > 1) {
    std::size_t chunk_size = (n + n_threads - 1) / n_threads;
    std::vector<std::thread> threads;
    for (std::size_t t = 0; t < n_threads; ++t) {
      std::size_t start_row = t * chunk_size;
      std::size_t end_row = std::min(start_row + chunk_size, n);
      threads.emplace_back(mmds_grad, std::cref(R_vec), std::cref(D_vec),
                           std::cref(Y_vec), eps, std::ref(gradient_vec),
                           start_row, end_row, n, d);
    }
    for (auto &thread : threads) {
      thread.join();
    }
  } else {
    mmds_grad(R_vec, D_vec, Y_vec, eps, gradient_vec, 0, n, n, d);
  }

  NumericMatrix gradient(n, d);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < d; ++j) {
      gradient(i, j) = gradient_vec[i * d + j];
    }
  }

  return gradient;
}
