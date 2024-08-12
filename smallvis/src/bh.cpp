#include "bh.h"
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix bh_tsne_gradient_cpp(IntegerVector indices, IntegerVector indptr,
                                   NumericVector P_data,
                                   NumericMatrix embedding, double theta = 0.5,
                                   double eps = 1e-16, int n_threads = 1) {
  std::vector<std::size_t> indices_cpp = as<std::vector<std::size_t>>(indices);
  std::vector<std::size_t> indptr_cpp = as<std::vector<std::size_t>>(indptr);
  std::vector<double> P_data_cpp = as<std::vector<double>>(P_data);

  std::vector<double> embedding_vec =
      as<std::vector<double>>(transpose(embedding));

  std::vector<double> gradient_cpp =
      smallvis::bh_tsne_gradient(indices_cpp, indptr_cpp, P_data_cpp,
                                 embedding_vec, theta, eps, n_threads);

  std::size_t n_samples = embedding.nrow();
  NumericMatrix gradient(2, n_samples, gradient_cpp.begin());
  return transpose(gradient);
}

// [[Rcpp::export]]
double bh_plogq_cpp(IntegerVector indices, IntegerVector indptr,
                    NumericVector P_data, NumericMatrix embedding,
                    double theta = 0.5, double eps = 1e-16, int n_threads = 1) {
  std::vector<std::size_t> indices_cpp = as<std::vector<std::size_t>>(indices);
  std::vector<std::size_t> indptr_cpp = as<std::vector<std::size_t>>(indptr);
  std::vector<double> P_data_cpp = as<std::vector<double>>(P_data);

  std::vector<double> embedding_vec =
      as<std::vector<double>>(transpose(embedding));

  return smallvis::bh_plogq(indices_cpp, indptr_cpp, P_data_cpp, embedding_vec,
                            theta, eps, n_threads);
}
