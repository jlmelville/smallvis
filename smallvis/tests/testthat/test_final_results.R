test_that("final costs match returned coordinates before the first epoch", {
  X <- matrix(sin(seq_len(96)^2), 24, 4)
  for (iter0_cost in c(FALSE, TRUE)) {
    res <- smallvis(
      X,
      perplexity = 4,
      max_iter = 10,
      epoch = 5,
      opt = list("l-bfgs", max_gr = 1),
      iter0_cost = iter0_cost,
      ret_extra = c("P", "Q"),
      epoch_callback = FALSE,
      verbose = FALSE
    )
    expect_equal(res$iter, 1)
    expect_length(res$itercosts, as.integer(iter0_cost))
    expect_equal(
      sum(res$costs),
      tsne_reference_cost(res$P, res$Y),
      tolerance = 1e-10
    )
    W <- 1 / (1 + as.matrix(stats::dist(res$Y))^2)
    diag(W) <- 0
    expect_equal(res$Q, W / sum(W), ignore_attr = TRUE, tolerance = 1e-10)
  }
})

test_that("final scoring does not reuse the previous epoch's cached point", {
  step <- opt_step
  local_mocked_bindings(opt_step = function(opt, cost_fn, Y, iter) {
    res <- step(opt, cost_fn, Y, iter)
    if (iter == 3) {
      res$opt$is_terminated <- TRUE
      res$opt$terminate <- list(what = "test stop after an epoch")
    }
    res
  })
  epoch_Y <- NULL
  res <- smallvis(
    iris10,
    perplexity = 4,
    Y_init = iris10_Y,
    max_iter = 8,
    epoch = 2,
    eta = 10,
    ret_extra = c("P", "Q"),
    verbose = FALSE,
    epoch_callback = function(Y, iter) if (iter == 2) epoch_Y <<- Y
  )
  expect_equal(res$iter, 3)
  expect_equal(names(res$itercosts), "2")
  expect_equal(unname(res$itercosts), tsne_reference_cost(res$P, epoch_Y))
  expect_gt(
    abs(tsne_reference_cost(res$P, res$Y) - unname(res$itercosts)),
    1e-6
  )
  expect_equal(
    sum(res$costs),
    tsne_reference_cost(res$P, res$Y),
    tolerance = 1e-10
  )
})

test_that("final scoring removes exaggeration without rewriting epoch history", {
  res <- smallvis(
    iris10,
    perplexity = 4,
    Y_init = iris10_Y,
    max_iter = 2,
    epoch = 2,
    exaggeration_factor = 12,
    stop_lying_iter = 10,
    ret_extra = "P",
    epoch_callback = FALSE,
    verbose = FALSE
  )
  expect_equal(sum(res$P), 1)
  expect_equal(
    sum(res$costs),
    tsne_reference_cost(res$P, res$Y),
    tolerance = 1e-10
  )
  expect_gt(abs(unname(res$itercosts) - sum(res$costs)), 1)
})

test_that("zero-iteration repeats and immediate L-BFGS stops have usable scores", {
  for (max_iter in c(0, 10)) {
    args <- list(
      X = iris10,
      perplexity = 4,
      nrep = 2,
      max_iter = max_iter,
      opt = list("l-bfgs", max_gr = 1),
      epoch = 5,
      epoch_callback = FALSE,
      verbose = FALSE
    )
    set.seed(42)
    all <- do.call(
      smallvis_rep,
      c(args, list(keep_all = TRUE, ret_extra = "P"))
    )
    scores <- vapply(
      all,
      function(res) tsne_reference_cost(res$P, res$Y),
      numeric(1)
    )
    expect_equal(all[[1]]$all_costs, scores, tolerance = 1e-10)
    set.seed(42)
    best <- do.call(smallvis_rep, c(args, list(ret_extra = FALSE)))
    expect_equal(best, all[[which.min(scores)]]$Y)
    set.seed(42)
    matrices <- do.call(
      smallvis_rep,
      c(args, list(keep_all = TRUE, ret_extra = FALSE))
    )
    expect_equal(matrices, lapply(all, `[[`, "Y"))
  }
})
