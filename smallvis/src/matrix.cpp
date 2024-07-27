#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <thread>
#include <vector>

using namespace Rcpp;

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
