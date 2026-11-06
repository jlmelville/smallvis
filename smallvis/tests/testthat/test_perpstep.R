test_that("an empty stepping schedule is equivalent to a direct call", {
  X <- matrix(sin(seq_len(96)^2), 24, 4)
  args <- list(
    X = X,
    perplexity = 10,
    max_iter = 4,
    ret_extra = "P",
    opt = list("adam", eta = 0.1),
    epoch_callback = FALSE,
    verbose = TRUE
  )
  set.seed(42)
  invisible(capture.output(
    stepped <- do.call(smallvis_perpstep, args),
    type = "message"
  ))
  set.seed(42)
  invisible(capture.output(direct <- do.call(smallvis, args), type = "message"))
  expect_equal(stepped$Y, direct$Y)
  expect_equal(stepped$P, direct$P)
  expect_equal(stepped$itercosts, direct$itercosts)
  expect_equal(stepped$costs, direct$costs)
})

test_that("preliminary settings do not leak into final default or explicit controls", {
  X <- matrix(sin(seq_len(96)^2), 24, 4)
  for (explicit in c(FALSE, TRUE)) {
    args <- list(
      X = X,
      perplexity = 3,
      max_iter = 12,
      ret_extra = "P",
      opt = list("adam", eta = 0.1),
      epoch_callback = FALSE,
      verbose = FALSE
    )
    if (explicit) {
      args$tol <- 1e-2
      args$epoch <- 2
      args$exaggeration_factor <- 3
      args$stop_lying_iter <- 2
    }
    set.seed(42)
    stepped <- do.call(smallvis_perpstep, c(args, list(step_iter = 4)))
    preliminary <- args
    preliminary$ret_extra <- FALSE
    preliminary$max_iter <- 2
    preliminary$tol <- 0
    preliminary$epoch <- 2
    preliminary$exaggeration_factor <- 1
    preliminary$late_exaggeration_factor <- 1
    set.seed(42)
    for (perp in c(8, 4)) {
      preliminary$perplexity <- perp
      preliminary$Y_init <- do.call(smallvis, preliminary)
    }
    args$Y_init <- preliminary$Y_init
    args$max_iter <- 8
    direct <- do.call(smallvis, args)
    expect_equal(stepped$Y, direct$Y)
    expect_equal(stepped$P, direct$P)
    expect_equal(stepped$itercosts, direct$itercosts)
    expect_equal(stepped$costs, direct$costs)
  }
})

test_that("stepping resolves and temporarily overrides nested method perplexity", {
  seen <- numeric()
  initialize <- sne_init
  local_mocked_bindings(sne_init = function(cost, X, perplexity, ...) {
    seen <<- c(seen, perplexity)
    initialize(cost, X, perplexity, ...)
  })
  X <- matrix(sin(seq_len(96)^2), 24, 4)
  res <- smallvis_perpstep(
    step_iter = 4,
    X = X,
    method = list("tsne", perplexity = 3),
    perplexity = 2,
    max_iter = 8,
    opt = list("adam", eta = 0.1),
    ret_extra = TRUE,
    epoch_callback = FALSE,
    verbose = FALSE
  )
  expect_equal(seen, c(8, 4, 3))
  expect_equal(res$perplexity, 3)
})

test_that("zero stepping iterations preserve a direct positional-input call", {
  X <- matrix(sin(seq_len(96)^2), 24, 4)
  args <- list(
    X,
    perplexity = 3,
    max_iter = 4,
    ret_extra = "P",
    opt = list("adam", eta = 0.1),
    epoch_callback = FALSE,
    verbose = FALSE
  )
  set.seed(42)
  stepped <- do.call(smallvis_perpstep, c(args, list(step_iter = 0)))
  set.seed(42)
  direct <- do.call(smallvis, args)
  expect_equal(stepped$Y, direct$Y)
  expect_equal(stepped$costs, direct$costs)
})
