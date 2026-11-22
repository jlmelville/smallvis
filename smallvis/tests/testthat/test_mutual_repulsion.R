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
  for (interval in list(0, -1, 1.5, Inf, NA_integer_, TRUE, c(1, 4), 2^31)) {
    expect_error(
      run(list("tsne", repulsion = "mutual", rebuild_interval = interval)),
      "rebuild_interval must be a positive integer"
    )
  }
  expect_error(
    run(list("tee", repulsion = "mutual", rebuild_interval = 4)),
    "Unknown method option"
  )
  expect_error(
    run(
      list("tsne", repulsion = "mutual", rebuild_interval = 4),
      opt = list("bfgs"),
      exaggeration_factor = 1
    ),
    "requires a matching objective and gradient"
  )
})

test_that("mutual refits update geometry and rebuild on the scheduled gradient", {
  i <- seq_len(512)
  y <- cbind((i %% 7) * 5 + sin(i), cos(i / 5))
  moved <- y[order(sin(i)), ]
  for (workers in c(1L, 2L)) {
    state <- mutual_workspace_cpp(4L)
    gradient <- function(y, theta = 0.35) {
      mutual_cauchy_cpp(y, theta, FALSE, workers, state)
    }
    fresh <- function(y, theta = 0.35) {
      mutual_cauchy_cpp(y, theta, FALSE, workers)
    }
    first <- gradient(y)
    expect_identical(first, fresh(y))
    # A rigid rotation preserves the cell openings, and rotates the force.
    rotated <- cbind(-y[, 2], y[, 1])
    second <- gradient(rotated)
    expect_equal(second$Z, first$Z, tolerance = 1e-12)
    expect_equal(
      second$force,
      cbind(-first$force[, 2], first$force[, 1]),
      tolerance = 1e-12
    )
    expect_equal(gradient(moved, 0), fresh(moved, 0), tolerance = 1e-12)
    expect_identical(
      mutual_cauchy_cpp(moved, 0.35, TRUE, workers, state)$Z,
      fresh(moved)$Z
    )
    fourth <- gradient(moved)
    expect_gt(max(abs(fourth$force - fresh(moved)$force)), 1e-8)
    expect_identical(gradient(moved), fresh(moved))
  }
})

test_that("each mutual t-SNE start owns its gradient tree", {
  i <- seq_len(128)
  args <- list(
    X = cbind(sin(i), cos(i), i),
    method = list("tsne", repulsion = "mutual", rebuild_interval = 4),
    perplexity = 4,
    nn = "exact",
    n_threads = 1L,
    eta = 10,
    max_iter = 7,
    epoch = 2,
    ret_extra = TRUE,
    verbose = FALSE
  )
  starts <- list()
  set.seed(42)
  runs <- do.call(
    smallvis_rep,
    c(
      args,
      list(
        nrep = 2,
        keep_all = TRUE,
        epoch_callback = function(Y, iter) {
          if (iter == 0) starts[[length(starts) + 1L]] <<- Y
        }
      )
    )
  )
  for (i in seq_along(runs)) {
    fit <- do.call(
      smallvis,
      c(args, list(Y_init = starts[[i]], epoch_callback = FALSE))
    )
    expect_identical(runs[[i]]$Y, fit$Y)
    expect_identical(runs[[i]]$costs, fit$costs)
  }
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
