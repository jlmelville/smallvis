test_that("Gaussian FFT normalization excludes the interpolated self mass", {
  Y <- unname(as.matrix(expand.grid(x = 3 * (0:4), y = 3 * (0:4)))[1:24, ])
  d2 <- as.matrix(stats::dist(Y))^2
  W <- exp(-d2)
  diag(W) <- 0
  log_mass <- rowSums(-d2 * W)
  sums <- fft_sne_sums(Y, 0.1, stats::fft, alpha = 0, log_weighted = TRUE)
  expect_equal(sums$Z, sum(W), tolerance = 0.005)
  expect_equal(sums$log_mass, unname(log_mass), tolerance = 0.005)
  expect_equal(
    sums$force,
    unname(Y * rowSums(W) - W %*% Y),
    tolerance = 0.005
  )

  expect_null(fft_sne_sums(10 * Y, 3, stats::fft, alpha = 0))
})

test_that("FFT SNE methods match exact BH on the same sparse affinities", {
  X <- cbind(seq_len(24), sin(seq_len(24)), cos(seq_len(24)))
  neighbors <- get_nn(stats::dist(X), k = 10L)
  Y <- cbind(sin(seq_len(24)), cos(seq_len(24)))
  Y[2L, ] <- Y[1L, ]
  separated <- as.matrix(expand.grid(x = 30 * (0:4), y = 30 * (0:4)))[1:24, ]
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
  cases <- list(
    list("hssne", alpha = 0.5),
    list("hssne", alpha = 1),
    list("dhssne", alpha = 0.5),
    list("ssne"),
    list("snerv", lambda = 0),
    list("snerv", lambda = 0.9),
    list("snerv", lambda = 1)
  )
  for (method in cases) {
    bh <- prepare(method, bh = TRUE, theta = 0)
    fft <- prepare(c(method, list(repulsion = "fft", grid_interval = 3)))
    expect_identical(fft$P, bh$P)
    coords <- list(Y)
    if (method[[1L]] %in% c("ssne", "snerv")) {
      coords <- c(coords, list(separated))
    }
    for (Z in coords) {
      expect_equal(
        cost_point(fft, Z)$pcost,
        cost_point(bh, Z)$pcost,
        tolerance = 1e-5
      )
      expect_equal(cost_grad(fft, Z)$G, cost_grad(bh, Z)$G, tolerance = 1e-5)
    }
  }
})

test_that("FFT dynamic HSSNE makes the same alpha choices as exact BH", {
  X <- cbind(seq_len(24), sin(seq_len(24)), cos(seq_len(24)))
  neighbors <- get_nn(stats::dist(X), k = 10L)
  Y <- cbind(sin(seq_len(24)), cos(seq_len(24))) * 3
  Y[2L, ] <- Y[1L, ]
  trace_alpha <- function(method, ...) {
    cost <- smallvis_run(
      X,
      method = method,
      neighbors = neighbors,
      perplexity = 3,
      scale = FALSE,
      prepare_only = TRUE,
      epoch_callback = FALSE,
      verbose = FALSE,
      ...
    )$cost_fn
    vapply(
      seq_len(4L),
      function(iter) {
        score <- cost_point(cost, Y)$pcost
        cost <<- cost$epoch(NULL, cost, iter, Y, score)$cost
        cost$alpha
      },
      numeric(1L)
    )
  }

  bh <- trace_alpha("dhssne", bh = TRUE, theta = 0)
  fft <- trace_alpha(list("dhssne", repulsion = "fft"))
  expect_true(any(bh != 0.5))
  expect_identical(fft, bh)
})
