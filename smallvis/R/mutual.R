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
  theta = 0.35,
  rebuild_interval = 1L
) {
  if (
    !is.numeric(rebuild_interval) ||
      length(rebuild_interval) != 1L ||
      !is.finite(rebuild_interval) ||
      rebuild_interval < 1 ||
      rebuild_interval > .Machine$integer.max ||
      rebuild_interval != trunc(rebuild_interval)
  ) {
    stop("rebuild_interval must be a positive integer")
  }
  cost <- cauchy_tsne_cost(
    perplexity,
    inp_kernel,
    nn,
    eps,
    n_threads,
    mutual_settings(theta, n_threads)
  )
  if (rebuild_interval > 1L) {
    cost$gradient_only <- TRUE
    cost$gr <- function(cost, Y) {
      # Prepared input can be shared by multiple starts. Each fit gets its
      # own tree on its first gradient, after input preparation.
      if (is.null(cost$mutual_workspace)) {
        cost$mutual_workspace <- mutual_workspace_cpp(rebuild_interval)
      }
      sums <- mutual_cauchy_cpp(
        Y,
        theta,
        FALSE,
        n_threads,
        cost$mutual_workspace
      )
      cauchy_tsne_gradient(cost, Y, sums, eps, n_threads)
    }
  }
  cost
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
