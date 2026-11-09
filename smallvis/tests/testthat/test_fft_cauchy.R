test_that("FFT Cauchy methods reuse saved affinities and agree with exact BH", {
  X <- cbind(seq_len(24), sin(seq_len(24)), cos(seq_len(24)))
  neighbors <- get_nn(stats::dist(X), k = 10L)
  Y <- cbind(sin(seq_len(24)), cos(seq_len(24)))
  Y[2L, ] <- Y[1L, ]
  local_mocked_bindings(get_nn = function(...) stop("neighbor search ran"))
  prepare <- function(method, ...) {
    smallvis_run(
      X,
      method = method,
      neighbors = neighbors,
      perplexity = 3,
      scale = FALSE,
      n_threads = 2,
      prepare_only = TRUE,
      epoch_callback = FALSE,
      verbose = FALSE,
      ...
    )$cost_fn
  }

  for (method in c("tee", "jssne", "hlsne", "chsne")) {
    bh <- prepare(method, bh = TRUE, theta = 0)
    fft <- prepare(list(method, repulsion = "fft"))
    expect_identical(fft$P, bh$P)
    for (factor in c(1, 4)) {
      fft_scaled <- start_exaggerating(fft, factor)
      bh_scaled <- start_exaggerating(bh, factor)
      for (coords in list(Y, matrix(0, nrow(Y), 2L))) {
        expect_equal(
          cost_grad(fft_scaled, coords)$G,
          cost_grad(bh_scaled, coords)$G,
          tolerance = 1e-5
        )
        expect_equal(
          cost_point(fft_scaled, coords)$pcost,
          cost_point(bh_scaled, coords)$pcost,
          tolerance = 1e-5
        )
      }
    }
    fit <- smallvis(
      X,
      method = list(method, repulsion = "fft"),
      neighbors = neighbors,
      perplexity = 3,
      scale = FALSE,
      Y_init = Y,
      eta = 1,
      max_iter = 2,
      ret_extra = "P",
      epoch_callback = FALSE,
      verbose = FALSE
    )
    expect_identical(fit$P, bh$P)
    expect_true(all(is.finite(fit$Y)) && is.finite(fit$costs))
  }
})

test_that("new FFT methods validate repulsion controls before graph setup", {
  local_mocked_bindings(get_nn = function(...) stop("neighbor search ran"))
  for (method in c("tee", "jssne", "hlsne", "chsne")) {
    fft <- list(method, repulsion = "fft")
    expect_error(smallvis(iris10, method = fft, bh = TRUE), "conflicts")
    expect_error(
      smallvis(iris10, method = fft, theta = 0),
      "theta applies only"
    )
    expect_error(
      smallvis(iris10, method = fft, opt = list("bfgs")),
      "requires a matching objective and gradient"
    )
  }
})
