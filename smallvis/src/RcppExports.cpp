// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// bh_tsne_gradient_cpp
NumericMatrix bh_tsne_gradient_cpp(IntegerVector indices, IntegerVector indptr, NumericVector P_data, NumericMatrix embedding, double theta, double eps, int n_threads);
RcppExport SEXP _smallvis_bh_tsne_gradient_cpp(SEXP indicesSEXP, SEXP indptrSEXP, SEXP P_dataSEXP, SEXP embeddingSEXP, SEXP thetaSEXP, SEXP epsSEXP, SEXP n_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type indices(indicesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type indptr(indptrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type P_data(P_dataSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type embedding(embeddingSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type n_threads(n_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(bh_tsne_gradient_cpp(indices, indptr, P_data, embedding, theta, eps, n_threads));
    return rcpp_result_gen;
END_RCPP
}
// bh_plogq_cpp
double bh_plogq_cpp(IntegerVector indices, IntegerVector indptr, NumericVector P_data, NumericMatrix embedding, double theta, double eps, int n_threads);
RcppExport SEXP _smallvis_bh_plogq_cpp(SEXP indicesSEXP, SEXP indptrSEXP, SEXP P_dataSEXP, SEXP embeddingSEXP, SEXP thetaSEXP, SEXP epsSEXP, SEXP n_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type indices(indicesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type indptr(indptrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type P_data(P_dataSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type embedding(embeddingSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type n_threads(n_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(bh_plogq_cpp(indices, indptr, P_data, embedding, theta, eps, n_threads));
    return rcpp_result_gen;
END_RCPP
}
// dist2_cpp
NumericMatrix dist2_cpp(NumericMatrix input, std::size_t n_threads);
RcppExport SEXP _smallvis_dist2_cpp(SEXP inputSEXP, SEXP n_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type input(inputSEXP);
    Rcpp::traits::input_parameter< std::size_t >::type n_threads(n_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(dist2_cpp(input, n_threads));
    return rcpp_result_gen;
END_RCPP
}
// dist_cpp
NumericMatrix dist_cpp(NumericMatrix input, std::size_t n_threads);
RcppExport SEXP _smallvis_dist_cpp(SEXP inputSEXP, SEXP n_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type input(inputSEXP);
    Rcpp::traits::input_parameter< std::size_t >::type n_threads(n_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(dist_cpp(input, n_threads));
    return rcpp_result_gen;
END_RCPP
}
// tweight_cpp
NumericMatrix tweight_cpp(NumericMatrix input, std::size_t n_threads);
RcppExport SEXP _smallvis_tweight_cpp(SEXP inputSEXP, SEXP n_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type input(inputSEXP);
    Rcpp::traits::input_parameter< std::size_t >::type n_threads(n_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(tweight_cpp(input, n_threads));
    return rcpp_result_gen;
END_RCPP
}
// d2_to_tweight_cpp
NumericMatrix d2_to_tweight_cpp(NumericMatrix dist_matrix, std::size_t n_threads);
RcppExport SEXP _smallvis_d2_to_tweight_cpp(SEXP dist_matrixSEXP, SEXP n_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type dist_matrix(dist_matrixSEXP);
    Rcpp::traits::input_parameter< std::size_t >::type n_threads(n_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(d2_to_tweight_cpp(dist_matrix, n_threads));
    return rcpp_result_gen;
END_RCPP
}
// tsne_grad_cpp
NumericMatrix tsne_grad_cpp(const NumericMatrix& P, const NumericMatrix& W, double Z, const NumericMatrix& Y, std::size_t n_threads);
RcppExport SEXP _smallvis_tsne_grad_cpp(SEXP PSEXP, SEXP WSEXP, SEXP ZSEXP, SEXP YSEXP, SEXP n_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type P(PSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type W(WSEXP);
    Rcpp::traits::input_parameter< double >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< std::size_t >::type n_threads(n_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(tsne_grad_cpp(P, W, Z, Y, n_threads));
    return rcpp_result_gen;
END_RCPP
}
// mmds_grad_cpp
NumericMatrix mmds_grad_cpp(const NumericMatrix& R, const NumericMatrix& D, const NumericMatrix& Y, double eps, std::size_t n_threads);
RcppExport SEXP _smallvis_mmds_grad_cpp(SEXP RSEXP, SEXP DSEXP, SEXP YSEXP, SEXP epsSEXP, SEXP n_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type R(RSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type D(DSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< std::size_t >::type n_threads(n_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(mmds_grad_cpp(R, D, Y, eps, n_threads));
    return rcpp_result_gen;
END_RCPP
}
// find_beta_knn_cpp
List find_beta_knn_cpp(const NumericMatrix& knn_distances, const IntegerMatrix& knn_indices, double perplexity, double tol, int max_tries, bool ret_sparse, std::size_t n_threads);
RcppExport SEXP _smallvis_find_beta_knn_cpp(SEXP knn_distancesSEXP, SEXP knn_indicesSEXP, SEXP perplexitySEXP, SEXP tolSEXP, SEXP max_triesSEXP, SEXP ret_sparseSEXP, SEXP n_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type knn_distances(knn_distancesSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type knn_indices(knn_indicesSEXP);
    Rcpp::traits::input_parameter< double >::type perplexity(perplexitySEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type max_tries(max_triesSEXP);
    Rcpp::traits::input_parameter< bool >::type ret_sparse(ret_sparseSEXP);
    Rcpp::traits::input_parameter< std::size_t >::type n_threads(n_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(find_beta_knn_cpp(knn_distances, knn_indices, perplexity, tol, max_tries, ret_sparse, n_threads));
    return rcpp_result_gen;
END_RCPP
}
// find_beta_cpp
List find_beta_cpp(const NumericMatrix& X, double perplexity, double tol, int max_tries, std::size_t n_threads);
RcppExport SEXP _smallvis_find_beta_cpp(SEXP XSEXP, SEXP perplexitySEXP, SEXP tolSEXP, SEXP max_triesSEXP, SEXP n_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type perplexity(perplexitySEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type max_tries(max_triesSEXP);
    Rcpp::traits::input_parameter< std::size_t >::type n_threads(n_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(find_beta_cpp(X, perplexity, tol, max_tries, n_threads));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_smallvis_bh_tsne_gradient_cpp", (DL_FUNC) &_smallvis_bh_tsne_gradient_cpp, 7},
    {"_smallvis_bh_plogq_cpp", (DL_FUNC) &_smallvis_bh_plogq_cpp, 7},
    {"_smallvis_dist2_cpp", (DL_FUNC) &_smallvis_dist2_cpp, 2},
    {"_smallvis_dist_cpp", (DL_FUNC) &_smallvis_dist_cpp, 2},
    {"_smallvis_tweight_cpp", (DL_FUNC) &_smallvis_tweight_cpp, 2},
    {"_smallvis_d2_to_tweight_cpp", (DL_FUNC) &_smallvis_d2_to_tweight_cpp, 2},
    {"_smallvis_tsne_grad_cpp", (DL_FUNC) &_smallvis_tsne_grad_cpp, 5},
    {"_smallvis_mmds_grad_cpp", (DL_FUNC) &_smallvis_mmds_grad_cpp, 5},
    {"_smallvis_find_beta_knn_cpp", (DL_FUNC) &_smallvis_find_beta_knn_cpp, 7},
    {"_smallvis_find_beta_cpp", (DL_FUNC) &_smallvis_find_beta_cpp, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_smallvis(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
