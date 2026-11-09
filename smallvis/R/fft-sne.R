fft_hssne <- function(
  perplexity,
  alpha = 0.5,
  inp_kernel = "nngaussian",
  symmetrize = "symmetric",
  eps = .Machine$double.eps,
  n_threads = 0,
  nn = "approximate",
  fft_backend = "stats",
  grid_interval = 1
) {
  alpha <- validate_hssne_alpha(alpha, "FFT")
  fft_normalized_sne(
    perplexity,
    inp_kernel,
    symmetrize,
    eps,
    n_threads,
    nn,
    alpha,
    fft_settings(fft_backend, grid_interval, n_threads)
  )
}

fft_dhssne <- function(
  perplexity,
  alpha = 0.5,
  inp_kernel = "nngaussian",
  symmetrize = "symmetric",
  eps = .Machine$double.eps,
  n_threads = 0,
  nn = "approximate",
  fft_backend = "stats",
  grid_interval = 1
) {
  cost <- fft_hssne(
    perplexity,
    alpha,
    inp_kernel,
    symmetrize,
    eps,
    n_threads,
    nn,
    fft_backend,
    grid_interval
  )
  cost$epoch <- dhssne_epoch
  cost
}

fft_ssne <- function(
  perplexity,
  inp_kernel = "nngaussian",
  symmetrize = "symmetric",
  n_threads = 0,
  nn = "approximate",
  fft_backend = "stats",
  grid_interval = 1
) {
  fft_normalized_sne(
    perplexity,
    inp_kernel,
    symmetrize,
    eps = 0,
    n_threads = n_threads,
    nn = nn,
    alpha = 0,
    fft = fft_settings(fft_backend, grid_interval, n_threads)
  )
}

fft_normalized_sne <- function(
  perplexity,
  inp_kernel,
  symmetrize,
  eps,
  n_threads,
  nn,
  alpha,
  fft
) {
  evaluate <- function(P, Y, alpha, gradient) {
    sums <- fft$sums(Y, mass_only = !gradient, alpha = alpha)
    if (is.null(sums)) {
      fn <- if (gradient) bh_sne_gradient_cpp else bh_sne_score_cpp
      return(fn(P@i, P@p, P@x, Y, 0, 0, n_threads, alpha, TRUE, 0))
    }
    if (gradient) {
      attraction <- sne_sparse_attraction_cpp(
        P@i,
        P@p,
        P@x,
        Y,
        n_threads,
        alpha
      )
      return(attraction - sums$force / sums$Z)
    }
    -sne_edge_loss(P, Y, log(sums$Z), alpha)
  }
  cost <- sparse_sne_cost(
    perplexity,
    inp_kernel,
    symmetrize,
    eps,
    n_threads,
    nn,
    alpha,
    normalize = TRUE,
    lambda = 0,
    score = function(P, Y, alpha) evaluate(P, Y, alpha, FALSE),
    gradient = function(P, Y, alpha) evaluate(P, Y, alpha, TRUE)
  )
  use_fft_cost(cost, fft)
}

fft_rklsne <- function(
  perplexity,
  inp_kernel = "nngaussian",
  symmetrize = "symmetric",
  eps = .Machine$double.eps,
  n_threads = 0,
  nn = "approximate",
  fft_backend = "stats",
  grid_interval = 1
) {
  fft_mixed_kl(
    perplexity,
    inp_kernel,
    symmetrize,
    eps,
    n_threads,
    nn,
    lambda = 0,
    gaussian = FALSE,
    fft = fft_settings(fft_backend, grid_interval, n_threads)
  )
}

fft_snerv <- function(
  perplexity,
  lambda = 0.9,
  inp_kernel = "nngaussian",
  symmetrize = "symmetric",
  eps = .Machine$double.eps,
  n_threads = 0,
  nn = "approximate",
  fft_backend = "stats",
  grid_interval = 1
) {
  fft_mixed_kl(
    perplexity,
    inp_kernel,
    symmetrize,
    eps,
    n_threads,
    nn,
    lambda,
    gaussian = TRUE,
    fft = fft_settings(fft_backend, grid_interval, n_threads)
  )
}

fft_mixed_kl <- function(
  perplexity,
  inp_kernel,
  symmetrize,
  eps,
  n_threads,
  nn,
  lambda,
  gaussian,
  fft
) {
  alpha <- if (gaussian) 0 else 1
  evaluate <- function(cost, Y, gradient) {
    P <- cost$P
    sums <- fft$sums(
      Y,
      mass_only = !gradient,
      alpha = alpha,
      log_weighted = TRUE
    )
    if (is.null(sums)) {
      return(bh_mixed_kl_cpp(
        P@i,
        P@p,
        cost$log_ratio,
        cost$p_excess,
        Y,
        0,
        n_threads,
        eps,
        lambda,
        cost$p_mass,
        cost$plogp,
        gaussian,
        gradient
      ))
    }
    force <- if (gradient) sums$force else NULL
    log_force <- if (gradient) sums$log_force else NULL
    fft_mixed_kl_cpp(
      P@i,
      P@p,
      cost$log_ratio,
      cost$p_excess,
      Y,
      sums$Z,
      sums$log_mass,
      force,
      log_force,
      n_threads,
      eps,
      lambda,
      cost$p_mass,
      cost$plogp,
      gaussian
    )
  }
  cost <- sparse_mixed_kl_cost(
    perplexity,
    inp_kernel,
    symmetrize,
    eps,
    n_threads,
    nn,
    lambda,
    evaluate,
    "FFT"
  )
  use_fft_cost(cost, fft)
}

sne_edge_loss <- function(P, Y, offset = 0, alpha = 1) {
  col <- rep(seq_len(ncol(P)), diff(P@p))
  row <- P@i + 1L
  dx <- Y[col, 1L] - Y[row, 1L]
  dy <- Y[col, 2L] - Y[row, 2L]
  d2 <- dx^2 + dy^2
  loss <- if (alpha == 0) d2 else log1p(alpha * d2) / alpha
  sum(P@x * (loss + offset))
}
