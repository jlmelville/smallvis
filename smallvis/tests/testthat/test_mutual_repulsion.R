test_that("mutual tree uses the public sparse input and exact limit for five objectives", {
  X <- cbind(seq_len(60), sin(seq_len(60)), cos(seq_len(60)))
  nn <- get_nn(stats::dist(X), k = 10L)
  Y <- cbind(sin(seq_len(60)), cos(seq_len(60)))
  Y[2, ] <- Y[1, ]
  run <- function(method, ...) {
    smallvis(
      X,
      method = method,
      neighbors = nn,
      perplexity = 3,
      scale = FALSE,
      Y_init = Y,
      eta = 0.5,
      max_iter = 2,
      ret_extra = "P",
      epoch_callback = FALSE,
      verbose = FALSE,
      ...
    )
  }
  for (method in c("tsne", "tee", "jssne", "hlsne", "chsne")) {
    options <- switch(
      method,
      tsne = list(eps = 0.1),
      tee = list(lambda = 0.07, symmetrize = "mutual"),
      jssne = list(eps = 1e-8),
      chsne = list(symmetrize = "fuzzy"),
      list()
    )
    bh <- run(
      c(list(method, repulsion = "bh"), options),
      theta = 0,
      n_threads = 2
    )
    mutual <- run(
      c(list(method, repulsion = "mutual"), options),
      theta = 0,
      n_threads = 2
    )
    expect_identical(mutual$P, bh$P)
    expect_equal(mutual$Y, bh$Y, tolerance = 1e-11)
    expect_equal(mutual$costs, bh$costs, tolerance = 1e-11)
    default <- run(list(method, repulsion = "mutual"))
    explicit <- run(list(method, repulsion = "mutual", theta = 0.35))
    expect_identical(default$Y, explicit$Y)
    expect_identical(default$costs, explicit$costs)
  }
})

test_that("hard cell expansions include both cells and exclude self", {
  # Well-separated cells exercise the Taylor expansion and its downward pass.
  i <- seq_len(128)
  y <- cbind((i %% 4) * 10 + 0.2 * sin(i), cos(i))
  ref <- mutual_cauchy_cpp(y, 0, FALSE)
  approx <- mutual_cauchy_cpp(y, 0.35, FALSE)
  expect_lt(abs(approx$Z / ref$Z - 1), 0.005)
  error <- sqrt(sum((ref$force - approx$force)^2) / sum(ref$force^2))
  expect_gt(error, 1e-10)
  expect_lt(error, 0.02)
  expect_equal(colSums(approx$force), c(0, 0), tolerance = 1e-12)
  expect_identical(mutual_cauchy_cpp(y, 0.35, TRUE)$Z, approx$Z)
  coincident <- mutual_cauchy_cpp(matrix(0, 64, 2), 0.35, FALSE)
  expect_identical(coincident$Z, 64 * 63)
  expect_true(all(coincident$force == 0))
})

test_that("mutual tree rejects incompatible controls before graph setup", {
  run <- function(method = list("tsne", repulsion = "mutual"), ...) {
    smallvis(
      iris10,
      method = method,
      epoch_callback = FALSE,
      verbose = FALSE,
      ...
    )
  }
  expect_error(run(bh = TRUE), "conflicts")
  expect_error(run(k = 3), "Mutual tree requires k = 2")
  expect_error(run(theta = -1), "finite and nonnegative")
  expect_error(
    run(list("umap", repulsion = "mutual")),
    "Mutual tree repulsion supports"
  )
  expect_error(
    run(list("tsne", repulsion = "mutual", grid_interval = 1)),
    "Unknown method option"
  )
})

test_that("parallel mutual pairs reproduce the serial sum and conserve force", {
  i <- seq_len(2048)
  y <- cbind((i %% 7) * 5 + sin(i), cos(i / 5))
  y[2L, ] <- y[1L, ]
  for (theta in c(0, 0.35)) {
    serial <- mutual_cauchy_cpp(y, theta, FALSE, 1L)
    parallel <- mutual_cauchy_cpp(y, theta, FALSE, 4L)
    expect_equal(parallel, serial, tolerance = 1e-11)
    expect_identical(parallel, mutual_cauchy_cpp(y, theta, FALSE, 4L))
    expect_identical(parallel$Z, mutual_cauchy_cpp(y, theta, TRUE, 4L)$Z)
    expect_lt(max(abs(colSums(parallel$force))), 1e-8)
  }
})
