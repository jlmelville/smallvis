fft_tee <- function(
  perplexity,
  inp_kernel = "nngaussian",
  symmetrize = "symmetric",
  lambda = 0.01,
  eps = .Machine$double.eps,
  n_threads = 0,
  nn = "approximate",
  fft_backend = "stats",
  grid_interval = 1,
  fft_reuse = NULL
) {
  validate_tee_lambda(lambda)
  fft <- fft_settings(fft_backend, grid_interval, n_threads, fft_reuse)
  cost <- sparse_sne_cost(
    perplexity,
    inp_kernel,
    symmetrize,
    eps,
    n_threads,
    nn,
    alpha = 1,
    normalize = FALSE,
    lambda = lambda,
    score = function(P, Y, alpha) {
      Z <- fft$sums(Y, mass_only = TRUE)$Z
      -sne_edge_loss(P, Y) - lambda * Z
    },
    gradient = function(P, Y, alpha) {
      sums <- fft$sums(Y)
      sne_sparse_attraction_cpp(P@i, P@p, P@x, Y, n_threads) -
        lambda * sums$force
    }
  )
  use_fft_cost(cost, fft)
}

fft_jssne <- function(
  perplexity,
  inp_kernel = "nngaussian",
  symmetrize = "symmetric",
  eps = .Machine$double.eps,
  n_threads = 0,
  nn = "approximate",
  fft_backend = "stats",
  grid_interval = 1,
  fft_reuse = NULL
) {
  sparse_divergence_cost(
    "jssne",
    perplexity,
    inp_kernel,
    symmetrize,
    n_threads,
    nn,
    eps = eps,
    fft = fft_settings(fft_backend, grid_interval, n_threads, fft_reuse)
  )
}

fft_hlsne <- function(
  perplexity,
  inp_kernel = "nngaussian",
  n_threads = 0,
  nn = "approximate",
  fft_backend = "stats",
  grid_interval = 1,
  fft_reuse = NULL
) {
  sparse_divergence_cost(
    "hlsne",
    perplexity,
    inp_kernel,
    "symmetric",
    n_threads,
    nn,
    fft = fft_settings(fft_backend, grid_interval, n_threads, fft_reuse)
  )
}

fft_chsne <- function(
  perplexity,
  inp_kernel = "nngaussian",
  symmetrize = "symmetric",
  n_threads = 0,
  nn = "approximate",
  fft_backend = "stats",
  grid_interval = 1,
  fft_reuse = NULL
) {
  sparse_divergence_cost(
    "chsne",
    perplexity,
    inp_kernel,
    symmetrize,
    n_threads,
    nn,
    fft = fft_settings(fft_backend, grid_interval, n_threads, fft_reuse)
  )
}
