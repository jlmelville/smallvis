test_that("AdaMax retains an elementwise recurrence and safeguards zero gradients", {
  # Arbitrary gradient sequences isolate the optimizer arithmetic from a cost.
  eta <- 0.3
  beta1 <- 0.8
  beta2 <- 0.7
  eps <- 0.01
  opt <- opt_init(adamax(eta, beta1, beta2, eps), 2, 2)
  m <- u <- matrix(0, 2, 2)
  gradients <- list(matrix(c(1, 0, 4, -2), 2), matrix(c(-3, 0, 2, -1), 2))
  for (step in seq_along(gradients)) {
    G <- gradients[[step]]
    expected <- matrix(0, 2, 2)
    for (i in 1:2) {
      for (j in 1:2) {
        m[i, j] <- beta1 * m[i, j] + (1 - beta1) * G[i, j]
        u[i, j] <- max(beta2 * u[i, j], abs(G[i, j]))
        expected[i, j] <- -eta * m[i, j] / ((1 - beta1^step) * (u[i, j] + eps))
      }
    }
    opt <- opt$upd(opt, NULL, matrix(0, 2, 2), step, G)
    expect_equal(opt$u, u)
    expect_equal(opt$uY, expected)
    expect_true(all(is.finite(opt$uY)))
  }
})

test_that("line searches reject unmatched forces before input calibration", {
  local_mocked_bindings(sne_init = function(...) stop("input calibration ran"))
  for (method in c("umap", "tumap", "largevis", "usne", "pstsne", "dhssne")) {
    expect_error(
      smallvis(
        iris10,
        method = method,
        perplexity = 4,
        opt = list("l-bfgs"),
        max_iter = 2,
        epoch_callback = FALSE,
        verbose = FALSE
      ),
      "requires a matching objective and gradient.*gradient-only optimizer"
    )
  }
})

test_that("unstabilized USNE matches its independent scalar derivative", {
  set.seed(42)
  X <- matrix(rnorm(24 * 4), 24)
  Y <- matrix(rnorm(48), 24)
  direction <- Y / sqrt(sum(Y * Y))
  for (min_dist in c(0.001, -0.5)) {
    cost <- cost_init(
      usne(perplexity = 4, gr_eps = 0, min_dist = min_dist),
      X,
      max_iter = 2
    )
    scalar <- function(Y) {
      W <- 1 / (1 + cost$a * as.matrix(stats::dist(Y))^(2 * cost$b))
      diag(W) <- 0
      Q <- W / sum(W)
      used <- cost$P > 0
      sum(cost$P[used] * log(cost$P[used] / Q[used]))
    }
    G <- cost_grad(cost_clear(cost), Y)$G
    expect_true(all(is.finite(G)))
    expect_equal(sum(cost_point(cost_clear(cost), Y)$pcost), scalar(Y))
    for (h in c(1e-3, 1e-4, 1e-5, 1e-6)) {
      numeric <- (scalar(Y + h * direction) - scalar(Y - h * direction)) /
        (2 * h)
      expect_equal(sum(G * direction), numeric, tolerance = 1e-5)
    }
    fit <- smallvis(
      X,
      method = list("usne", gr_eps = 0, min_dist = min_dist),
      scale = FALSE,
      perplexity = 4,
      Y_init = Y,
      opt = list("l-bfgs"),
      max_iter = 2,
      ret_extra = TRUE,
      epoch_callback = FALSE,
      verbose = FALSE
    )
    expect_true(all(is.finite(fit$Y)))
    expect_true(all(is.finite(fit$costs)))
  }
})

test_that("USNE handles coincident pairs only where its kernel is differentiable", {
  X <- matrix(sin(seq_len(96)^2), 24, 4)
  Y <- cbind(sin(1:24), cos(1:24))
  Y[1:2, ] <- 0
  cost <- cost_init(usne(perplexity = 4, gr_eps = 0), X, max_iter = 2)
  expect_gt(unname(cost$b), 0.5)
  G <- cost_grad(cost_clear(cost), Y)$G
  expect_true(all(is.finite(G)))
  direction <- matrix(sin(seq_len(48)), 24)
  direction <- direction / sqrt(sum(direction^2))
  scalar <- function(Y) sum(cost_point(cost_clear(cost), Y)$pcost)
  for (h in c(1e-4, 1e-5, 1e-6)) {
    numeric <- (scalar(Y + h * direction) - scalar(Y - h * direction)) /
      (2 * h)
    expect_equal(sum(G * direction), numeric, tolerance = 1e-5)
  }
  run <- function(gr_eps, min_dist = 0.001) {
    smallvis(
      X,
      method = list("usne", gr_eps = gr_eps, min_dist = min_dist),
      scale = FALSE,
      perplexity = 4,
      Y_init = Y,
      opt = list("adam", eta = 0.001),
      max_iter = 2,
      epoch_callback = FALSE,
      verbose = FALSE
    )
  }
  expect_true(all(is.finite(run(0))))
  expect_error(run(0, -0.5), "not differentiable at coincident coordinates")
  expect_true(all(is.finite(run(0.1, -0.5))))
})

test_that("USNE and pseudo-separable t-SNE retain gradient-only optimization", {
  X <- matrix(sin(seq_len(96)^2), 24, 4)
  Y <- cbind(sin(1:24), cos(1:24))
  for (method in c("usne", "pstsne")) {
    fit <- smallvis(
      X,
      method = method,
      scale = FALSE,
      perplexity = 4,
      Y_init = Y,
      opt = list("adam", eta = 0.001),
      max_iter = 3,
      epoch = 1,
      tol = 0,
      min_cost = -Inf,
      ret_extra = TRUE,
      epoch_callback = FALSE,
      verbose = FALSE
    )
    expect_equal(fit$iter, 3)
    expect_true(all(is.finite(fit$Y)))
    expect_true(all(is.finite(fit$costs)))
  }
})

test_that("unstabilized forces match directional derivatives at several step sizes", {
  X <- matrix(sin(seq_len(96)^2), 24, 4)
  Y <- cbind(sin(1:24), cos(1:24))
  direction <- cbind(cos(2 * (1:24)), sin(3 * (1:24)))
  direction <- direction / sqrt(sum(direction^2))
  for (method in c("umap", "tumap", "largevis")) {
    cost <- do.call(get(method), list(perplexity = 4, gr_eps = 0))
    cost <- cost_init(cost, X, max_iter = 2)
    analytic <- sum(cost_grad(cost_clear(cost), Y)$G * direction)
    scalar <- function(Y) sum(cost_point(cost_clear(cost), Y)$pcost)
    for (h in c(1e-4, 1e-5, 1e-6)) {
      numeric <- (scalar(Y + h * direction) - scalar(Y - h * direction)) /
        (2 * h)
      expect_equal(analytic, numeric, tolerance = 1e-5, info = method)
    }
    res <- smallvis(
      X,
      method = list(method, gr_eps = 0),
      Y_init = Y,
      perplexity = 4,
      opt = list("l-bfgs"),
      max_iter = 2,
      epoch_callback = FALSE,
      verbose = FALSE
    )
    expect_true(all(is.finite(res)))
    stabilized <- smallvis(
      X,
      method = method,
      Y_init = Y,
      perplexity = 4,
      max_iter = 1,
      epoch_callback = FALSE,
      verbose = FALSE
    )
    expect_true(all(is.finite(stabilized)))
  }
  # This normalized variant replaces the stabilized force it inherits.
  expect_no_error(smallvis(
    X,
    method = "ntumap",
    perplexity = 4,
    Y_init = Y,
    opt = list("l-bfgs"),
    max_iter = 1,
    epoch_callback = FALSE,
    verbose = FALSE
  ))
})

test_that("t-EE exaggeration scales the attractive force and restores P", {
  X <- matrix(sin(seq_len(96)^2), 24, 4)
  Y <- cbind(sin(1:24), cos(1:24))
  lambda <- 0.1
  cost <- cost_init(tee(perplexity = 4, lambda = lambda), X, max_iter = 3)
  P <- cost$P
  D2 <- as.matrix(stats::dist(Y))^2
  W <- 1 / (1 + D2)
  diag(W) <- 0
  for (factor in c(1, 12)) {
    exaggerated <- start_exaggerating(cost, factor)
    actual <- cost_grad(cost_clear(exaggerated), Y)$G
    expected <- matrix(0, 24, 2)
    for (i in 1:24) {
      for (j in 1:24) {
        expected[i, ] <- expected[i, ] +
          4 /
            sum(P) *
            W[i, j] *
            (factor * P[i, j] - lambda * W[i, j]) *
            (Y[i, ] - Y[j, ])
      }
    }
    expect_equal(actual, expected, tolerance = 1e-10)
    restored <- stop_exaggerating(exaggerated, factor)
    expect_equal(restored$P, P)
    expect_equal(
      cost_grad(cost_clear(restored), Y)$G,
      cost_grad(cost_clear(cost), Y)$G
    )
  }
  run <- function(factor) {
    smallvis(
      X,
      method = "tee",
      perplexity = 4,
      Y_init = Y,
      max_iter = 2,
      exaggeration_factor = factor,
      stop_lying_iter = 2,
      opt = list("steepd", eta = 0.1),
      epoch_callback = FALSE,
      verbose = FALSE
    )
  }
  expect_gt(max(abs(run(12) - run(1))), 1e-3)
})
