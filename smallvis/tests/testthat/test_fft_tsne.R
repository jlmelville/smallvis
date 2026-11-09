test_that("explicit exact and Barnes-Hut selectors preserve the legacy routes", {
  run <- function(...) {
    smallvis(
      iris10,
      perplexity = 2,
      Y_init = "spca",
      max_iter = 2,
      ret_extra = "P",
      epoch_callback = FALSE,
      verbose = FALSE,
      ...
    )
  }
  exact <- run()
  explicit <- run(method = list("tsne", repulsion = "exact"))
  expect_identical(explicit$P, exact$P)
  expect_identical(explicit$Y, exact$Y)

  bh <- run(bh = TRUE, theta = 0, nn = "exact")
  explicit <- run(
    method = list("tsne", repulsion = "bh", theta = 0),
    bh = TRUE,
    nn = "exact"
  )
  expect_identical(explicit$P, bh$P)
  expect_identical(explicit$Y, bh$Y)
})

test_that("incompatible FFT settings fail before searching for neighbors", {
  local_mocked_bindings(get_nn = function(...) stop("neighbor search ran"))
  fft <- list("tsne", repulsion = "fft")
  run <- function(method = fft, ...) {
    smallvis(
      iris10,
      method = method,
      epoch_callback = FALSE,
      verbose = FALSE,
      ...
    )
  }
  expect_error(run(bh = TRUE), "conflicts with repulsion")
  expect_error(run(bh = NA), "bh must be TRUE or FALSE")
  expect_error(run(list("tsne", repulsion = "exact"), bh = TRUE), "conflicts")
  expect_error(run(theta = 0), "theta applies only")
  expect_error(run(c(fft, list(theta = 0.5))), "theta applies only")
  expect_error(run(k = 3), "FFT repulsion requires k = 2")
  expect_error(run(inp_kernel = "gaussian"), "nngaussian.*knn")
  expect_error(
    run(opt = list("bfgs")),
    "requires a matching objective and gradient"
  )
  expect_error(run(c(fft, list(grid_interval = 0))), "finite positive")
  expect_error(run(c(fft, list(fft_backend = "missing"))), "fft_backend must")
  expect_error(
    run(list("tsne", fft_backend = "stats")),
    "Unknown method option"
  )
  expect_error(run(list("tsne", repulsion = "auto")), "repulsion must")
  expect_error(run(list("asne", repulsion = "fft")), "Unknown method option")
  expect_error(run(neighbors = list(), nn = "exact"), "Leave nn unset")
})

test_that("FFT t-SNE uses saved sparse affinities and approximates exact updates", {
  X <- cbind(seq_len(24), sin(seq_len(24)), cos(seq_len(24)))
  neighbors <- get_nn(stats::dist(X), k = 10L)
  Y <- cbind(sin(seq_len(24)), cos(seq_len(24)))
  local_mocked_bindings(get_nn = function(...) stop("neighbor search ran"))
  run <- function(method, ...) {
    smallvis(
      X,
      method = method,
      perplexity = 3,
      neighbors = neighbors,
      scale = FALSE,
      Y_init = Y,
      eta = 1,
      max_iter = 1,
      ret_extra = "P",
      epoch_callback = FALSE,
      verbose = FALSE,
      ...
    )
  }
  exact <- run(list("tsne", repulsion = "bh", theta = 0))
  fft <- run(list("tsne", repulsion = "fft"))
  expect_s4_class(fft$P, "dgCMatrix")
  expect_identical(fft$P, exact$P)
  expect_equal(fft$Y, exact$Y, tolerance = 1e-5)
  expect_equal(fft$costs, tsne_reference_cost(fft$P, fft$Y), tolerance = 1e-5)

  if (requireNamespace("fftwtools", quietly = TRUE)) {
    fftw <- run(
      list("tsne", repulsion = "fft", fft_backend = "fftwtools"),
      n_threads = 2L
    )
    expect_equal(fftw$Y, fft$Y, tolerance = 1e-12)
    expect_equal(fftw$costs, fft$costs, tolerance = 1e-12)
  }

  if (requireNamespace("zufft", quietly = TRUE)) {
    real_fft <- run(
      list("tsne", repulsion = "fft", fft_backend = "zufft"),
      n_threads = 2L
    )
    expect_equal(real_fft$Y, fft$Y, tolerance = 1e-12)
    expect_equal(real_fft$costs, fft$costs, tolerance = 1e-12)
  }
})

test_that("FFT padding keeps interpolation nodes fixed on an odd-sized grid", {
  Y <- rbind(c(0, 0), c(1, 1), c(202, 202))
  grid <- fft_sne_grid(Y, 1)
  expect_equal(grid$width, 1)
  expect_equal(grid$side, 1215)
  expect_equal(dim(grid$S), c((3 * 202)^2, nrow(Y)))

  skip_if_not_installed("zufft")
  baseline <- fft_settings("stats", 1, 0)$sums(Y)
  real_fft <- fft_settings("zufft", 1, 0)$sums(Y)
  expect_equal(real_fft$Z, baseline$Z, tolerance = 1e-12)
  # Cancellation makes the distant point's tiny force sensitive to roundoff.
  force_scale <- max(abs(baseline$force))
  expect_lt(max(abs(real_fft$force - baseline$force)) / force_scale, 1e-11)
})

test_that("FFT t-SNE evaluates coincident coordinates and guards wide grids", {
  run <- function(Y) {
    smallvis(
      iris10,
      method = list("tsne", repulsion = "fft"),
      nn = "exact",
      perplexity = 2,
      Y_init = Y,
      max_iter = 0,
      ret_extra = "P",
      epoch_callback = FALSE,
      verbose = FALSE
    )
  }
  Y <- cbind(rep(3, nrow(iris10)), rep(5, nrow(iris10)))
  result <- run(Y)
  expect_equal(
    result$costs,
    tsne_reference_cost(result$P, Y),
    tolerance = 1e-12
  )
  Y[1L, 1L] <- 10000
  # The final score is evaluated after recentering even at max_iter = 0.
  expect_error(run(Y), "padded grid.*Increase grid_interval")
})
