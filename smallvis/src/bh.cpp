#include "bh.h"
#include <Rcpp.h>

using namespace Rcpp;

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
      smallvis::bh_tsne_gradient(indices_cpp, indptr_cpp, P_data_cpp,
                                 embedding_vec, theta, eps, n_threads);

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

  return smallvis::bh_plogq(indices_cpp, indptr_cpp, P_data_cpp, embedding_vec,
                            theta, eps, n_threads);
}
