# Piecewise interpolation of the Cauchy kernel and its square, as in
# Linderman et al. (2019), doi:10.1038/s41592-018-0308-4.
fft_tsne <- function(
  perplexity,
  inp_kernel = "nngaussian",
  nn = "approximate",
  eps = .Machine$double.eps,
  n_threads = 0,
  fft_backend = "stats",
  grid_interval = 1,
  fft_reuse = NULL
) {
  if (
    !is.character(inp_kernel) ||
      length(inp_kernel) != 1L ||
      !inp_kernel %in% c("nngaussian", "knn")
  ) {
    stop("FFT t-SNE requires inp_kernel = 'nngaussian' or 'knn'")
  }
  nn <- match.arg(nn, c("approximate", "exact"))
  if (
    !is.numeric(n_threads) ||
      length(n_threads) != 1L ||
      !is.finite(n_threads) ||
      n_threads < 0 ||
      n_threads > .Machine$integer.max ||
      n_threads != trunc(n_threads)
  ) {
    stop("n_threads must be a nonnegative integer")
  }
  if (
    !is.numeric(eps) ||
      is.complex(eps) ||
      length(eps) != 1L ||
      !is.finite(eps) ||
      eps < 0
  ) {
    stop("FFT t-SNE requires finite nonnegative eps")
  }
  fft <- fft_settings(fft_backend, grid_interval, n_threads, fft_reuse)
  cost <- sparse_tsne_cost(
    perplexity,
    inp_kernel,
    nn,
    n_threads,
    eps,
    pfn = function(cost, Y) {
      Z <- fft$sums(Y, mass_only = TRUE)$Z
      cost$pcost <- cost$plogp + sne_edge_loss(cost$P, Y, log(Z + eps))
      cost
    },
    gr = function(cost, Y) {
      repulsion <- fft$sums(Y)
      P <- cost$P
      attraction <- sne_sparse_attraction_cpp(P@i, P@p, P@x, Y, n_threads)
      cost$G <- 4 * (attraction - repulsion$force / (repulsion$Z + eps))
      cost
    }
  )
  use_fft_cost(cost, fft)
}
