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
  expect_error(run(c(fft, list(fft_reuse = NA))), "fft_reuse must")
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

  reuse <- run(list("tsne", repulsion = "fft", fft_reuse = TRUE))
  expect_identical(reuse$Y, fft$Y)
  expect_identical(reuse$costs, fft$costs)

  if (requireNamespace("fftwtools", quietly = TRUE)) {
    fftw <- run(
      list("tsne", repulsion = "fft", fft_backend = "fftwtools"),
      n_threads = 2L
    )
    expect_equal(fftw$Y, fft$Y, tolerance = 1e-12)
    expect_equal(fftw$costs, fft$costs, tolerance = 1e-12)
  }

  if (requireNamespace("zufft", quietly = TRUE)) {
    # Reuse gives identical coordinates, so check the settings at the native call.
    settings <- NULL
    native <- zufft_tsne_sums_cpp
    local_mocked_bindings(zufft_tsne_sums_cpp = function(
      workspace,
      points,
      interval,
      fast_boxes,
      workers,
      mass_only,
      reuse,
      hybrid
    ) {
      settings <<- list(workers, reuse)
      native(
        workspace,
        points,
        interval,
        fast_boxes,
        workers,
        mass_only,
        reuse,
        hybrid
      )
    })
    zufft <- list("tsne", repulsion = "fft", fft_backend = "zufft")
    one_worker <- run(zufft, n_threads = 0L)
    expect_identical(settings, list(1L, TRUE))
    expect_equal(one_worker$Y, fft$Y, tolerance = 1e-12)
    expect_equal(one_worker$costs, fft$costs, tolerance = 1e-12)

    two_workers <- run(zufft, n_threads = 2L)
    expect_identical(settings, list(2L, TRUE))
    expect_equal(two_workers$Y, one_worker$Y, tolerance = 1e-12)
    expect_equal(two_workers$costs, one_worker$costs, tolerance = 1e-12)

    separate <- run(c(zufft, list(fft_reuse = FALSE)), n_threads = 0L)
    expect_identical(settings, list(1L, FALSE))
    expect_identical(separate$Y, one_worker$Y)
    expect_identical(separate$costs, one_worker$costs)
  }
})

test_that("FFT keeps all nine stencil entries for the self-mass correction", {
  grid <- fft_sne_grid(rbind(c(0, 0), c(50, 2), c(0.5, 0.5)), 1)
  S <- grid$S
  expect_identical(diff(S@p), rep(9L, 3L))
  # The last point lands exactly on the center node of the first box.
  expect_identical(tail(S@i, 9), c(0:2, 150:152, 300:302))
  expect_equal(tail(S@x, 9), as.numeric(seq_len(9) == 5L))
  impulse <- matrix(0, grid$side[1L], grid$side[2L])
  impulse[1L, 1L] <- 1
  expect_equal(fft_self_mass(grid, impulse), Matrix::colSums(S^2))
})

test_that("FFT crops unused boxes and pads rectangular axes independently", {
  Y <- rbind(c(0, 50), c(1, 51), c(202, 67))
  grid <- fft_sne_grid(Y, 1)
  expect_equal(grid$width, 1)
  expect_equal(grid$side, c(1215, 108))
  expect_equal(dim(grid$S), c((3 * 202) * (3 * 18), nrow(Y)))

  baseline <- fft_settings("stats", 1, 0)$sums(Y)
  # Cancellation makes the distant point's tiny force sensitive to roundoff.
  force_scale <- max(abs(baseline$force))
  for (backend in c("fftwtools", "zufft")) {
    if (requireNamespace(backend, quietly = TRUE)) {
      result <- fft_settings(backend, 1, 0)$sums(Y)
      expect_equal(result$Z, baseline$Z, tolerance = 1e-12)
      expect_lt(max(abs(result$force - baseline$force)) / force_scale, 1e-11)
    }
  }
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
  Y[1L, ] <- 10000
  # The final score is evaluated after recentering even at max_iter = 0.
  expect_error(run(Y), "padded grid.*Increase grid_interval")
})

test_that("FFT spacing and t-SNE updates are invariant to translation", {
  X <- matrix(sin(seq_len(96)^2), 24)
  Y <- cbind(sin(1:24), cos(1:24))
  translated <- sweep(Y, 2L, c(1000, -300), "+")
  expect_equal(fft_sne_grid(translated, 1)$width, fft_sne_grid(Y, 1)$width)
  expect_equal(fft_sne_grid(translated, 1)$side, fft_sne_grid(Y, 1)$side)

  W2 <- (1 + as.matrix(stats::dist(Y))^2)^-2
  diag(W2) <- 0
  direct <- Y * rowSums(W2) - W2 %*% Y
  fft <- fft_settings("stats", 1, 1)
  for (Y_init in list(Y, translated)) {
    error <- max(abs(fft$sums(Y_init)$force - direct)) / max(abs(direct))
    expect_lt(error, 2e-5)
  }

  run <- function(Y_init, method) {
    smallvis(
      X,
      method = method,
      Y_init = Y_init,
      nn = "exact",
      perplexity = 4,
      scale = FALSE,
      opt = list("steepd", eta = 0.5),
      max_iter = 1,
      epoch_callback = FALSE,
      verbose = FALSE
    )
  }
  exact <- run(Y, list("tsne", repulsion = "bh", theta = 0))
  method <- list("tsne", repulsion = "fft")
  original <- run(Y, method)
  shifted <- run(translated, method)
  expect_lt(max(abs(shifted - exact)), 1e-5)
  expect_equal(shifted, original, tolerance = 1e-10)
})

test_that("fused zufft t-SNE refreshes its grid and preserves the Cauchy fields", {
  skip_if_not_installed("zufft")
  generic <- fft_settings("zufft", 0.75, 1)$sums
  Y <- cbind(sin(1:50), cos(1:50))
  wide <- sweep(Y, 2, c(40, 3), "*")
  flat <- cbind(wide[, 1], 42)
  for (reuse in c(FALSE, TRUE)) {
    fused <- zufft_tsne_sums(0.75, 2L, reuse)
    for (points in list(Y, wide, flat, Y)) {
      ref <- generic(points)
      observed <- fused(points)
      expect_equal(observed, ref, tolerance = 1e-11)
      expect_equal(fused(points, mass_only = TRUE), list(Z = observed$Z))
    }
    expect_equal(fused(matrix(2, 4, 2)), list(Z = 12, force = matrix(0, 4, 2)))
    expect_error(fused(Y * 1e6), "padded grid.*Increase grid_interval")
    expect_error(fused(Y * Inf), "coordinates must be finite")
  }
})

test_that("hybrid t-SNE routes gradient and score on either side of 14 cells per point", {
  skip_if_not_installed("zufft")
  # Both line layouts occupy a 300 by 6 padded grid. Adding one point changes
  # the ratio from just above 14 to just below it.
  line <- cbind(seq(0, 1, length.out = 129), 0)
  short <- line[-2, , drop = FALSE]
  wide <- line * 1e6
  tree <- mutual_settings(0.1, 2L)$sums
  fft <- zufft_tsne_sums(1, 2L, FALSE)
  hybrid <- zufft_tsne_sums(1, 2L, FALSE, tree)
  for (mass_only in c(FALSE, TRUE)) {
    expected <- fft(line, mass_only)
    expect_false(identical(expected, tree(line, mass_only)))
    expect_identical(hybrid(line, mass_only), expected)
    expect_identical(hybrid(short, mass_only), tree(short, mass_only))
    expect_identical(hybrid(wide, mass_only), tree(wide, mass_only))
    expect_identical(hybrid(line, mass_only), expected)
  }
  # A 1,050 by 6 grid fits exactly at 14 cells per point until zufft pads it.
  padded <- cbind(seq(0, 175, length.out = 450), 0)
  expect_identical(hybrid(padded), tree(padded))
  # 2,048^2 cells fit the per-point budget here, but exceed the absolute cap.
  large <- cbind(rep(c(0, 340), 160000), rep(c(0, 340), 160000))
  capped <- zufft_tsne_sums(1, 1L, TRUE, function(...) "tree")
  expect_identical(capped(large), "tree")
})

test_that("hybrid is an opt-in two-dimensional t-SNE method", {
  run <- function(method = list("tsne", repulsion = "hybrid"), ...) {
    smallvis(
      iris10,
      method = method,
      perplexity = 2,
      nn = "exact",
      Y_init = "spca",
      max_iter = 2,
      n_threads = 2L,
      ret_extra = "P",
      epoch_callback = FALSE,
      verbose = FALSE,
      ...
    )
  }
  # The small data takes the mutual route with or without optional zufft.
  hybrid <- run(theta = 0.1)
  pure <- run(list("tsne", repulsion = "mutual", theta = 0.1))
  expect_identical(hybrid$Y, pure$Y)
  expect_identical(hybrid$costs, pure$costs)
  expect_error(run(k = 3), "Hybrid repulsion requires k = 2")
  expect_error(
    run(opt = list("bfgs")),
    "requires a matching objective and gradient"
  )
  expect_error(run(bh = TRUE), "conflicts")
  expect_error(run(theta = -1), "finite and nonnegative")
  expect_error(run(list("tee", repulsion = "hybrid")), "only t-SNE")
  expect_error(
    run(list("tsne", repulsion = "hybrid", fft_reuse = NA)),
    "fft_reuse"
  )
})

test_that("hybrid uses the available evaluator for a compact grid", {
  y <- cbind(seq(0, 1, length.out = 150), 0)
  run <- function(method) {
    smallvis(
      y,
      method = method,
      perplexity = 2,
      nn = "exact",
      scale = FALSE,
      Y_init = y,
      n_threads = 2L,
      max_iter = 1,
      eta = 1,
      ret_extra = "P",
      epoch_callback = FALSE,
      verbose = FALSE
    )
  }
  hybrid <- run(list("tsne", repulsion = "hybrid"))
  method <- if (requireNamespace("zufft", quietly = TRUE)) {
    list("tsne", repulsion = "fft", fft_backend = "zufft")
  } else {
    list("tsne", repulsion = "mutual")
  }
  pure <- run(method)
  expect_identical(hybrid$Y, pure$Y)
  expect_identical(hybrid$costs, pure$costs)
})
