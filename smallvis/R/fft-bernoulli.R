fft_largevis <- function(
  perplexity,
  inp_kernel = "nngaussian",
  symmetrize = "symmetric",
  gamma = 1,
  gr_eps = 0.1,
  normalize = TRUE,
  eps = 1e-9,
  row_weight = NULL,
  n_threads = 0,
  nn = "approximate",
  fft_backend = "stats",
  grid_interval = 1,
  fft_reuse = NULL
) {
  fft <- fft_settings(fft_backend, grid_interval, n_threads, fft_reuse)
  cost <- bhlargevis(
    perplexity,
    inp_kernel,
    symmetrize,
    gamma,
    gr_eps,
    normalize,
    eps,
    row_weight,
    n_threads,
    theta = 0,
    nn = nn
  )
  fft_bernoulli_cost(cost, fft, "largevis", eps, gr_eps, gamma, n_threads)
}

fft_tumap <- function(
  perplexity,
  inp_kernel = "nngaussian",
  symmetrize = "umap",
  gr_eps = 0.1,
  eps = 1e-9,
  row_weight = NULL,
  n_threads = 0,
  nn = "approximate",
  fft_backend = "stats",
  grid_interval = 1,
  fft_reuse = NULL
) {
  fft <- fft_settings(fft_backend, grid_interval, n_threads, fft_reuse)
  cost <- bhtumap(
    perplexity,
    inp_kernel,
    symmetrize,
    gr_eps,
    eps,
    row_weight,
    n_threads,
    theta = 0,
    nn = nn
  )
  fft_bernoulli_cost(cost, fft, "tumap", eps, gr_eps, 1, n_threads)
}

fft_umap <- function(
  perplexity,
  inp_kernel = "skd",
  symmetrize = "umap",
  spread = 1,
  min_dist = 0.001,
  gr_eps = 0.1,
  eps = 1e-9,
  row_weight = NULL,
  n_threads = 0,
  nn = "approximate",
  fft_backend = "stats",
  grid_interval = 1,
  fft_reuse = NULL
) {
  fft <- fft_settings(fft_backend, grid_interval, n_threads, fft_reuse)
  cost <- bhumap(
    perplexity,
    inp_kernel,
    symmetrize,
    spread,
    min_dist,
    gr_eps,
    eps,
    row_weight,
    n_threads,
    theta = 0,
    nn = nn
  )
  fft_bernoulli_cost(cost, fft, "umap", eps, gr_eps, 1, n_threads)
}

fft_bernoulli_cost <- function(
  cost,
  fft,
  method,
  eps,
  gr_eps,
  gamma,
  n_threads
) {
  exact_score <- cost$pfn
  exact_gradient <- cost$gr
  evaluate <- function(cost, Y, gradient) {
    grid <- fft_sne_grid(Y, fft$grid_interval)
    if (is.null(grid)) {
      fn <- if (gradient) exact_gradient else exact_score
      return(fn(cost, Y))
    }
    a <- if (method == "umap") cost$a else 1
    b <- if (method == "umap") cost$b else 1
    # The log-complement objective is sharper than the stabilized force.
    radius <- (if (gradient) 1.5 else 2) * grid$width
    if (gradient) {
      # On a fine grid, the stabilized force is already smooth at grid scale.
      stabilization <- if (method == "largevis") gr_eps else gr_eps + eps
      scale <- min(sqrt(stabilization), a^(-0.5 / b))
      if (grid$width <= 0.1 * scale) {
        radius <- 0
      }
    }
    kernel <- fft_bernoulli_grid_cpp(
      grid$d2,
      method,
      eps,
      gr_eps,
      a,
      b,
      radius,
      gradient,
      n_threads
    )
    convolve <- fft_grid_convolver(grid, fft$transform, mass_only = !gradient)
    fields <- convolve(kernel, if (gradient) 1:3 else 1L)
    force <- if (gradient) fft_grid_force(grid, fields) else NULL
    mass <- if (gradient) 0 else sum(fields[, 1L] - fft_self_mass(grid, kernel))
    P <- cost$P
    result <- fft_bernoulli_finish_cpp(
      P@i,
      P@p,
      P@x,
      Y,
      force,
      mass,
      method,
      eps,
      gr_eps,
      a,
      b,
      gamma,
      radius,
      n_threads
    )
    if (gradient) {
      cost$G <- 4 * result
    } else {
      cost$pcost <- if (method == "largevis") -result else cost$Cp - result
    }
    cost
  }
  cost$pfn <- function(cost, Y) evaluate(cost, Y, gradient = FALSE)
  cost$gr <- function(cost, Y) evaluate(cost, Y, gradient = TRUE)
  use_fft_cost(cost, fft)
}
