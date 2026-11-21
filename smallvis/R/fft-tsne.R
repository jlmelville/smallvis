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
  engine <- fft_settings(fft_backend, grid_interval, n_threads, fft_reuse)
  if (fft_backend == "zufft") {
    engine$sums <- zufft_tsne_sums(grid_interval, n_threads, fft_reuse)
  }
  cauchy_tsne_cost(
    perplexity,
    inp_kernel,
    nn,
    eps,
    n_threads,
    engine
  )
}

# Each fit owns its buffers. Coordinates and kernel values change every update.
zufft_tsne_sums <- function(grid_interval, n_threads, fft_reuse) {
  workspace <- zufft_workspace_cpp()
  workers <- max(1L, n_threads)
  reuse <- is.null(fft_reuse) || isTRUE(fft_reuse)
  function(Y, mass_only = FALSE) {
    zufft_tsne_sums_cpp(
      workspace,
      Y,
      grid_interval,
      fft_fast_boxes,
      workers,
      mass_only,
      reuse
    )
  }
}

cauchy_tsne_cost <- function(
  perplexity,
  inp_kernel,
  nn,
  eps,
  n_threads,
  engine
) {
  label <- engine$label
  if (
    !is.character(inp_kernel) ||
      length(inp_kernel) != 1L ||
      !inp_kernel %in% c("nngaussian", "knn")
  ) {
    stop(label, " t-SNE requires inp_kernel = 'nngaussian' or 'knn'")
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
    stop(label, " t-SNE requires finite nonnegative eps")
  }
  cost <- sparse_tsne_cost(
    perplexity,
    inp_kernel,
    nn,
    n_threads,
    eps,
    pfn = function(cost, Y) {
      Z <- engine$sums(Y, mass_only = TRUE)$Z
      P <- cost$P
      cost$pcost <- cost$plogp +
        cauchy_edge_loss_cpp(P@i, P@p, P@x, Y, log(Z + eps), n_threads)
      cost
    },
    gr = function(cost, Y) {
      repulsion <- engine$sums(Y)
      P <- cost$P
      attraction <- sne_sparse_attraction_cpp(P@i, P@p, P@x, Y, n_threads)
      cost$G <- 4 * (attraction - repulsion$force / (repulsion$Z + eps))
      cost
    }
  )
  use_repulsion_cost(cost, engine)
}
