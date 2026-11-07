test_that("dynamic HSSNE never selects an alpha below its positive bound", {
  set.seed(42)
  X <- matrix(rnorm(24 * 4), 24)
  Y <- matrix(rnorm(48), 24) * 0.1
  seen <- numeric()
  fit <- smallvis(
    X,
    method = list("dhssne", alpha = 0.1),
    scale = FALSE,
    perplexity = 4,
    Y_init = Y,
    opt = list("adam", eta = 1e-5),
    max_iter = 6,
    epoch = 1,
    tol = 0,
    min_cost = -Inf,
    ret_extra = "alpha",
    epoch_callback = function(Y, iter, cost, cost_fn) {
      seen[iter + 1] <<- cost_fn$alpha
    },
    verbose = FALSE
  )
  expect_equal(fit$iter, 6)
  expect_length(seen, 7)
  expect_true(all(is.finite(seen)))
  expect_true(all(seen >= 1e-8))
  expect_gte(fit$alpha, 1e-8)
  expect_true(all(is.finite(fit$Y)))
  expect_true(all(is.finite(fit$costs)))
})

test_that("the first evaluation after an alpha change uses fresh weights", {
  set.seed(42)
  X <- matrix(rnorm(24 * 4), 24)
  Y <- matrix(rnorm(48), 24)
  old <- cost_init(dhssne(perplexity = 4, alpha = 0.5), X, max_iter = 2)
  old <- cost_point(cost_grad(old, Y), Y)
  changed <- do_epoch(opt_create(list("adam")), old, 1, Y, sum(old$pcost))$cost
  expect_equal(changed$alpha, 0.6)
  expect_equal(changed$apow, -1 / 0.6)
  expected <- cost_init(hssne(perplexity = 4, alpha = 0.6), X, max_iter = 2)
  expected <- cost_point(cost_grad(expected, Y), Y)
  # Derived data must be either invalidated or already correct for the new alpha.
  for (field in c("W", "Z", "G", "pcost")) {
    if (!is.null(changed[[field]])) {
      expect_equal(changed[[field]], expected[[field]], tolerance = 1e-12)
    }
  }
  # Deliberately avoid cost_clear(): the epoch transition owns invalidation.
  actual <- cost_point(cost_grad(changed, Y), Y)
  expect_equal(actual$G, expected$G, tolerance = 1e-12)
  expect_equal(actual$W, expected$W, tolerance = 1e-12)
  expect_equal(actual$Z, expected$Z, tolerance = 1e-12)
  expect_equal(actual$pcost, expected$pcost, tolerance = 1e-12)
  expect_equal(actual$P, old$P)
})

test_that("retaining alpha preserves the current objective and gradient", {
  set.seed(42)
  X <- matrix(rnorm(24 * 4), 24)
  Y <- matrix(rnorm(48), 24)
  current <- cost_init(dhssne(perplexity = 4, alpha = 0.5), X, max_iter = 2)
  # Construct a known optimum at the current alpha, independently of the cache.
  W <- (1 + 0.5 * as.matrix(stats::dist(Y))^2)^(-2)
  diag(W) <- 0
  current$P <- W / sum(W)
  current <- cost_cache_input(current)
  current <- cost_point(cost_grad(current, Y), Y)
  next_cost <- do_epoch(
    opt_create(list("adam")),
    current,
    1,
    Y,
    sum(current$pcost)
  )$cost
  expect_equal(next_cost$alpha, 0.5)
  actual <- cost_point(cost_grad(next_cost, Y), Y)
  expect_equal(actual$G, current$G, tolerance = 1e-12)
  expect_equal(actual$W, current$W, tolerance = 1e-12)
  expect_equal(sum(actual$pcost), 0, tolerance = 1e-12)
})

test_that("fixed-alpha HSSNE retains function-based optimization", {
  fit <- smallvis(
    matrix(sin(seq_len(96)^2), 24, 4),
    method = list("hssne", alpha = 0.5),
    perplexity = 4,
    Y_init = cbind(sin(1:24), cos(1:24)),
    opt = list("l-bfgs"),
    max_iter = 2,
    ret_extra = TRUE,
    epoch_callback = FALSE,
    verbose = FALSE
  )
  expect_true(all(is.finite(fit$Y)))
  expect_true(all(is.finite(fit$costs)))
})
