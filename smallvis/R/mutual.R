mutual_settings <- function(theta, n_threads) {
  if (
    !is.numeric(theta) ||
      is.complex(theta) ||
      length(theta) != 1L ||
      !is.finite(theta) ||
      theta < 0
  ) {
    stop("Mutual tree theta must be finite and nonnegative")
  }
  list(
    repulsion = "mutual",
    label = "Mutual tree",
    sums = function(Y, mass_only = FALSE) {
      mutual_cauchy_cpp(Y, theta, mass_only, n_threads)
    }
  )
}

mutual_tsne <- function(
  perplexity,
  inp_kernel = "nngaussian",
  nn = "approximate",
  eps = .Machine$double.eps,
  n_threads = 0,
  theta = 0.35
) {
  cauchy_tsne_cost(
    perplexity,
    inp_kernel,
    nn,
    eps,
    n_threads,
    mutual_settings(theta, n_threads)
  )
}

mutual_tee <- function(
  perplexity,
  inp_kernel = "nngaussian",
  symmetrize = "symmetric",
  lambda = 0.01,
  eps = .Machine$double.eps,
  n_threads = 0,
  nn = "approximate",
  theta = 0.35
) {
  cauchy_tee_cost(
    perplexity,
    inp_kernel,
    symmetrize,
    lambda,
    eps,
    n_threads,
    nn,
    mutual_settings(theta, n_threads)
  )
}

mutual_jssne <- function(
  perplexity,
  inp_kernel = "nngaussian",
  symmetrize = "symmetric",
  eps = .Machine$double.eps,
  n_threads = 0,
  nn = "approximate",
  theta = 0.35
) {
  sparse_divergence_cost(
    "jssne",
    perplexity,
    inp_kernel,
    symmetrize,
    n_threads,
    nn,
    eps = eps,
    engine = mutual_settings(theta, n_threads)
  )
}

mutual_hlsne <- function(
  perplexity,
  inp_kernel = "nngaussian",
  n_threads = 0,
  nn = "approximate",
  theta = 0.35
) {
  sparse_divergence_cost(
    "hlsne",
    perplexity,
    inp_kernel,
    "symmetric",
    n_threads,
    nn,
    engine = mutual_settings(theta, n_threads)
  )
}

mutual_chsne <- function(
  perplexity,
  inp_kernel = "nngaussian",
  symmetrize = "symmetric",
  n_threads = 0,
  nn = "approximate",
  theta = 0.35
) {
  sparse_divergence_cost(
    "chsne",
    perplexity,
    inp_kernel,
    symmetrize,
    n_threads,
    nn,
    engine = mutual_settings(theta, n_threads)
  )
}
