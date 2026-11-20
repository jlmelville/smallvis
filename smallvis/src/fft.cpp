#include <Rcpp.h>
#include <limits>

using namespace Rcpp;

// Each point is one sparse column. Visit its nine nodes in increasing row
// order so Matrix can use the indices and weights without sorting triplets.
// [[Rcpp::export(rng = false)]]
List fft_stencil_cpp(const NumericVector &first, int stride,
                     const NumericMatrix &wx, const NumericMatrix &wy) {
  const int n = wx.nrow();
  const int max_index = std::numeric_limits<int>::max();
  if (first.size() != n || wy.nrow() != n || wx.ncol() != 3 || wy.ncol() != 3 ||
      n > max_index / 9 || stride < 3 || stride > (max_index - 2) / 2) {
    stop("Invalid FFT stencil dimensions");
  }
  IntegerVector rows(no_init(9 * n));
  NumericVector weights(no_init(9 * n));
  for (int p = 0; p < n; ++p) {
    if (!(first[p] >= 0 && first[p] <= max_index - 2 * stride - 2)) {
      stop("Invalid FFT stencil offset");
    }
    const int start = static_cast<int>(first[p]);
    for (int y = 0; y < 3; ++y) {
      for (int x = 0; x < 3; ++x) {
        const int k = 9 * p + 3 * y + x;
        rows[k] = start + y * stride + x;
        weights[k] = wx(p, x) * wy(p, y);
      }
    }
  }
  return List::create(_["i"] = rows, _["x"] = weights);
}
