#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <thread>
#include <vector>

#include "threads.h"

using namespace Rcpp;

void mmds_grad(const std::vector<double> &R, const std::vector<double> &D,
               const std::vector<double> &Y, double eps,
               std::vector<double> &gradient, std::size_t start,
               std::size_t end, std::size_t n) {

  for (std::size_t i = start; i < end; ++i) {
    const std::size_t i_d = i * 2;
    const std::size_t i_n = i * n;
    for (std::size_t j = 0; j < n; ++j) {
      if (i == j) {
        continue;
      }
      const std::size_t ij = i_n + j;
      double k_ij = 4.0 * (D[ij] - R[ij]) / (D[ij] + eps);
      gradient[i_d] += k_ij * (Y[i_d] - Y[j * 2]);
      gradient[i_d + 1] += k_ij * (Y[i_d + 1] - Y[j * 2 + 1]);
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
             std::size_t start_row, std::size_t end_row, std::size_t n) {
  for (std::size_t i = start_row; i < end_row; ++i) {
    const std::size_t i_d = i * 2;
    const std::size_t i_n = i * n;
    for (std::size_t j = 0; j < n; ++j) {
      const std::size_t j_d = j * 2;
      double diff_x = data[i_d] - data[j_d];
      double diff_y = data[i_d + 1] - data[j_d + 1];
      dist_matrix[i_n + j] = 1.0 / (1.0 + diff_x * diff_x + diff_y * diff_y);
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
  // this can be used for the high dimensional data, do NOT assume d = 2!

  std::size_t n = input.nrow();
  std::size_t d = input.ncol();

  std::vector<double> transposed_data =
      as<std::vector<double>>(transpose(input));

  std::vector<double> dist_matrix(n * n, 0.0);
  auto worker = [&](std::size_t start, std::size_t end) {
    d2(transposed_data, dist_matrix, start, end, n, d);
  };
  smallvis::parallel_for(n, n_threads, worker);

  return NumericMatrix(n, n, dist_matrix.begin());
}

// [[Rcpp::export]]
NumericMatrix dist_cpp(NumericMatrix input, std::size_t n_threads = 1) {
  // this can be used for the high dimensional data, do NOT assume d = 2!

  std::size_t n = input.nrow();
  std::size_t d = input.ncol();

  std::vector<double> transposed_data =
      as<std::vector<double>>(transpose(input));

  std::vector<double> dist_matrix(n * n, 0.0);
  auto worker = [&](std::size_t start, std::size_t end) {
    dist(transposed_data, dist_matrix, start, end, n, d);
  };
  smallvis::parallel_for(n, n_threads, worker);

  return NumericMatrix(n, n, dist_matrix.begin());
}

// [[Rcpp::export]]
NumericMatrix tweight_cpp(NumericMatrix input, std::size_t n_threads = 1) {
  std::size_t n = input.nrow();

  std::vector<double> transposed_data =
      as<std::vector<double>>(transpose(input));

  std::vector<double> dist_matrix(n * n, 0.0);
  auto worker = [&](std::size_t start, std::size_t end) {
    tweight(transposed_data, dist_matrix, start, end, n);
  };
  smallvis::parallel_for(n, n_threads, worker);

  return NumericMatrix(n, n, dist_matrix.begin());
}

// [[Rcpp::export]]
NumericMatrix d2_to_tweight_cpp(NumericMatrix dist_matrix,
                                std::size_t n_threads) {
  std::size_t n = dist_matrix.nrow();

  std::vector<double> dist_matrix_vec(dist_matrix.begin(), dist_matrix.end());
  std::vector<double> transformed_matrix(n * n, 0.0);

  auto worker = [&](std::size_t start, std::size_t end) {
    d2_to_tweight(dist_matrix_vec, transformed_matrix, start, end);
  };
  smallvis::parallel_for(n, n_threads, worker);

  NumericMatrix result(n, n, transformed_matrix.begin());

  return result;
}

// [[Rcpp::export]]
NumericMatrix tsne_grad_cpp(const NumericMatrix &P, const NumericMatrix &W,
                            double Z, const NumericMatrix &Y,
                            std::size_t n_threads) {
  std::size_t n = Y.nrow();

  std::vector<double> P_vec(P.begin(), P.end());
  std::vector<double> W_vec(W.begin(), W.end());

  std::vector<double> Y_vec = as<std::vector<double>>(transpose(Y));

  std::vector<double> gradient_vec(n * 2, 0.0);
  auto worker = [&](std::size_t start, std::size_t end) {
    tsne_grad(P_vec, W_vec, Z, Y_vec, gradient_vec, start, end, n);
  };
  smallvis::parallel_for(n, n_threads, worker);

  NumericMatrix gradient(2, n, gradient_vec.begin());
  return transpose(gradient);
}

// [[Rcpp::export]]
NumericMatrix mmds_grad_cpp(const NumericMatrix &R, const NumericMatrix &D,
                            const NumericMatrix &Y, double eps,
                            std::size_t n_threads) {
  std::size_t n = Y.nrow();
  std::size_t d = Y.ncol();

  std::vector<double> R_vec(R.begin(), R.end());
  std::vector<double> D_vec(D.begin(), D.end());

  std::vector<double> Y_vec = as<std::vector<double>>(transpose(Y));

  std::vector<double> gradient_vec(n * d, 0.0);

  auto worker = [&](std::size_t start, std::size_t end) {
    mmds_grad(R_vec, D_vec, Y_vec, eps, gradient_vec, start, end, n);
  };
  smallvis::parallel_for(n, n_threads, worker);

  NumericMatrix gradient(2, n, gradient_vec.begin());
  return transpose(gradient);
}
