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
  cauchy_tee_cost(
    perplexity,
    inp_kernel,
    symmetrize,
    lambda,
    eps,
    n_threads,
    nn,
    fft_settings(fft_backend, grid_interval, n_threads, fft_reuse)
  )
}

cauchy_tee_cost <- function(
  perplexity,
  inp_kernel,
  symmetrize,
  lambda,
  eps,
  n_threads,
  nn,
  engine
) {
  validate_tee_lambda(lambda)
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
      Z <- engine$sums(Y, mass_only = TRUE)$Z
      -cauchy_edge_loss_cpp(P@i, P@p, P@x, Y, 0, n_threads) - lambda * Z
    },
    gradient = function(P, Y, alpha) {
      sums <- engine$sums(Y)
      sne_sparse_attraction_cpp(P@i, P@p, P@x, Y, n_threads) -
        lambda * sums$force
    }
  )
  use_repulsion_cost(cost, engine)
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
    engine = fft_settings(fft_backend, grid_interval, n_threads, fft_reuse)
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
    engine = fft_settings(fft_backend, grid_interval, n_threads, fft_reuse)
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
    engine = fft_settings(fft_backend, grid_interval, n_threads, fft_reuse)
  )
}
