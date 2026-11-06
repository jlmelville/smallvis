test_that("approximate-neighbor repeats share one input graph and rank actual scores", {
  X <- matrix(sin(seq_len(144)^2), 24, 6)
  calls <- 0L
  search <- get_nn
  local_mocked_bindings(get_nn = function(X, k, exact, n_threads, verbose) {
    calls <<- calls + 1L
    expect_false(exact)
    search(X, k, exact, n_threads, verbose)
  })
  set.seed(42)
  runs <- smallvis_rep(
    nrep = 3,
    keep_all = TRUE,
    X = X,
    perplexity = 4,
    bh = TRUE,
    theta = 0,
    nn = "approximate",
    n_threads = 1,
    opt = list("adam", eta = 0.05),
    max_iter = 4,
    ret_extra = "P",
    epoch_callback = FALSE,
    verbose = FALSE
  )
  expect_equal(calls, 1L)
  scores <- vapply(
    runs,
    function(res) tsne_reference_cost(res$P, res$Y),
    numeric(1)
  )
  for (res in runs) {
    expect_equal(res$P, runs[[1]]$P, tolerance = 0)
    expect_equal(
      sum(res$costs),
      tsne_reference_cost(res$P, res$Y),
      tolerance = 1e-10
    )
    expect_equal(res$all_costs, scores, tolerance = 1e-10)
    expect_equal(res$best_rep, which.min(scores))
  }
})

test_that("repeats prepare preprocessing once and reset optimizer and exaggeration state", {
  X <- matrix(sin(seq_len(144)^2), 24, 6)
  scale_calls <- 0L
  scale_original <- scale_input
  local_mocked_bindings(scale_input = function(...) {
    scale_calls <<- scale_calls + 1L
    scale_original(...)
  })
  starts <- list()
  args <- list(
    X = X,
    perplexity = 4,
    pca = TRUE,
    initial_dims = 3,
    opt = list("adam", eta = 0.1),
    max_iter = 4,
    exaggeration_factor = 3,
    stop_lying_iter = 1,
    late_exaggeration_factor = 2,
    start_late_lying_iter = 3,
    ret_extra = "P",
    verbose = FALSE
  )
  set.seed(42)
  runs <- do.call(
    smallvis_rep,
    c(
      args,
      list(
        nrep = 3,
        keep_all = TRUE,
        epoch_callback = function(Y, iter) {
          if (iter == 0) starts[[length(starts) + 1L]] <<- Y
        }
      )
    )
  )
  expect_equal(scale_calls, 1L)
  expect_length(starts, 3)
  expect_gt(max(abs(starts[[1]] - starts[[2]])), 1e-6)
  for (i in seq_along(runs)) {
    direct <- do.call(
      smallvis,
      c(args, list(Y_init = starts[[i]], epoch_callback = FALSE))
    )
    expect_equal(runs[[i]]$Y, direct$Y, tolerance = 1e-12)
    expect_equal(runs[[i]]$P, direct$P, tolerance = 1e-12)
    expect_equal(sum(runs[[i]]$P), 1)
  }
  scores <- vapply(
    runs,
    function(res) tsne_reference_cost(res$P, res$Y),
    numeric(1)
  )
  set.seed(42)
  best <- do.call(smallvis_rep, c(args, list(nrep = 3, epoch_callback = FALSE)))
  expect_equal(best$Y, runs[[which.min(scores)]]$Y)
  expect_equal(best$all_costs, scores, tolerance = 1e-10)
})

test_that("custom constructors retain per-run preparation", {
  calls <- 0L
  constructor <- function(perplexity) {
    calls <<- calls + 1L
    tsne(perplexity)
  }
  runs <- smallvis_rep(
    nrep = 2,
    keep_all = TRUE,
    X = iris10,
    method = constructor,
    perplexity = 4,
    max_iter = 0,
    ret_extra = TRUE,
    epoch_callback = FALSE,
    verbose = FALSE
  )
  expect_equal(calls, 2L)
  expect_identical(runs[[1]]$method, constructor)
  expect_length(runs[[1]]$all_costs, 2)
})
