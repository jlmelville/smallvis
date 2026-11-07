test_that("native dense gradients preserve caller-owned buffers", {
  # A private call makes ownership of P and W observable before R wrappers can
  # hide a native write through copy-on-modify.
  n <- 31
  Y <- cbind(sin(seq_len(n)), cos(seq_len(n)))
  W <- 1 / (1 + as.matrix(stats::dist(Y))^2)
  diag(W) <- 0
  P <- outer(seq_len(n), seq_len(n), function(i, j) 1 / (1 + abs(i - j)))
  diag(P) <- 0
  P <- P / sum(P)
  alias_P <- P
  alias_W <- W
  before <- unserialize(serialize(list(P = P, W = W, Y = Y), NULL))
  Z <- sum(W)
  expected <- matrix(0, n, 2)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      expected[i, ] <- expected[i, ] +
        4 * W[i, j] * (P[i, j] - W[i, j] / Z) * (Y[i, ] - Y[j, ])
    }
  }
  for (threads in c(1L, 2L)) {
    actual <- tsne_grad_cpp(P, W, Z, Y, n_threads = threads)
    expect_equal(actual, expected, tolerance = 1e-12)
    expect_identical(P, before$P)
    expect_identical(alias_P, before$P)
    expect_identical(W, before$W)
    expect_identical(alias_W, before$W)
    expect_identical(Y, before$Y)
  }
})

test_that("native MDS gradients borrow inputs without modifying aliases", {
  # Direct calls expose writes that R's copy-on-modify could hide in a wrapper.
  n <- 31
  X <- matrix(sin(seq_len(n * 4)^2), n, 4)
  Y <- cbind(sin(seq_len(n)), cos(seq_len(n)))
  Y[2, ] <- Y[1, ]
  R <- unname(as.matrix(stats::dist(X)))
  D <- unname(as.matrix(stats::dist(Y)))
  eps <- 0.01
  alias_R <- R
  alias_D <- D
  alias_Y <- Y
  before <- unserialize(serialize(list(R = R, D = D, Y = Y), NULL))
  expected <- matrix(0, n, 2)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      expected[i, ] <- expected[i, ] +
        4 * (D[i, j] - R[i, j]) / (D[i, j] + eps) * (Y[i, ] - Y[j, ])
    }
  }
  for (threads in c(1L, 2L)) {
    actual <- mmds_grad_cpp(R, D, Y, eps, n_threads = threads)
    expect_equal(actual, expected, tolerance = 1e-12)
    expect_identical(R, before$R)
    expect_identical(alias_R, before$R)
    expect_identical(D, before$D)
    expect_identical(alias_D, before$D)
    expect_identical(Y, before$Y)
    expect_identical(alias_Y, before$Y)
    # The same protected matrix may supply both read-only input buffers.
    expect_identical(
      mmds_grad_cpp(D, D, Y, eps, n_threads = threads),
      matrix(0, n, 2)
    )
  }
})

test_that("native Student weights own their result and preserve the diagonal", {
  # Check the raw result's lifetime and unit diagonal before t-SNE zeros it.
  Y <- cbind(c(0, 1, -2, 0, 4), c(0, -3, 0.5, 0, 2))
  alias_Y <- Y
  before <- unserialize(serialize(Y, NULL))
  expected <- 1 / (1 + unname(as.matrix(stats::dist(Y)))^2)
  for (threads in c(1L, 2L)) {
    actual <- tweight_cpp(Y, n_threads = threads)
    expect_type(actual, "double")
    expect_identical(dim(actual), c(5L, 5L))
    expect_equal(actual, expected, tolerance = 1e-15)
    expect_identical(diag(actual), rep(1, nrow(Y)))
    invisible(tweight_cpp(Y + 1, n_threads = threads))
    invisible(gc())
    expect_equal(actual, expected, tolerance = 1e-15)
    expect_identical(Y, before)
    expect_identical(alias_Y, before)

    cost <- cost_update(
      tsne(perplexity = 2, use_cpp = TRUE, n_threads = threads),
      Y
    )
    diag(expected) <- 0
    expect_equal(cost$W, expected, tolerance = 1e-15)
    expect_equal(cost$Z, sum(expected), tolerance = 1e-15)
    diag(expected) <- 1
  }
})

test_that("dense calibration preserves directed rows and input ownership", {
  n <- 47
  X <- matrix(sin(seq_len(n * 6)^2), n, 6)
  Y <- cbind(sin(seq_len(n)), cos(seq_len(n)))
  before <- unserialize(serialize(X, NULL))
  D2 <- as.matrix(stats::dist(X))^2
  reference <- NULL
  for (threads in c(1L, 2L)) {
    res <- smallvis(
      X,
      Y_init = Y,
      perplexity = 4,
      scale = FALSE,
      use_cpp = TRUE,
      n_threads = threads,
      max_iter = 0,
      ret_extra = c("V", "P", "beta"),
      epoch_callback = FALSE,
      verbose = FALSE
    )
    expected <- exp(-D2 * res$beta)
    diag(expected) <- 0
    expect_gt(max(abs(expected - t(expected))), 0.01)
    expect_equal(res$V, expected, tolerance = 1e-12, ignore_attr = TRUE)
    expect_equal(perpp(res$V), rep(4, n), tolerance = 1e-4)
    expect_identical(X, before)
    if (is.null(reference)) {
      reference <- res
    } else {
      expect_identical(res$V, reference$V)
      expect_identical(res$P, reference$P)
      expect_identical(res$beta, reference$beta)
    }
  }
})

test_that("reused calibration workspace handles interleaved failures and successes", {
  X <- matrix(sin(seq_len(32 * 4)^2), 32, 4)
  X[seq(1, 15, 2), ] <- 0
  one <- find_beta_cpp(X, perplexity = 4, n_threads = 1)
  two <- find_beta_cpp(X, perplexity = 4, n_threads = 2)
  expect_gt(one$bad_perp, 0)
  expect_lt(one$bad_perp, nrow(X))
  expect_identical(one, two)
  expect_equal(diag(one$W), rep(0, nrow(X)))
  expect_true(all(is.finite(one$W)))
  expect_true(all(rowSums(one$W) > 0))
})

test_that("complete dense native runs match independent fixed-affinity updates", {
  n <- 47
  X <- matrix(sin(seq_len(n * 6)^2), n, 6)
  Y <- cbind(sin(seq_len(n)), cos(seq_len(n)))
  before <- unserialize(serialize(list(X = X, Y = Y), NULL))
  initial <- smallvis(
    X,
    Y_init = Y,
    perplexity = 4,
    scale = FALSE,
    use_cpp = TRUE,
    max_iter = 0,
    ret_extra = "P",
    epoch_callback = FALSE,
    verbose = FALSE
  )
  expected <- Y
  for (step in 1:4) {
    W <- 1 / (1 + as.matrix(stats::dist(expected))^2)
    diag(W) <- 0
    Z <- sum(W)
    G <- matrix(0, n, 2)
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        G[i, ] <- G[i, ] +
          4 *
            W[i, j] *
            (initial$P[i, j] - W[i, j] / Z) *
            (expected[i, ] - expected[j, ])
      }
    }
    expected <- expected - 0.1 * G
    expected <- sweep(expected, 2, colMeans(expected))
  }
  for (threads in c(1L, 2L)) {
    actual <- smallvis(
      X,
      Y_init = Y,
      perplexity = 4,
      scale = FALSE,
      use_cpp = TRUE,
      n_threads = threads,
      max_iter = 4,
      opt = list("steepd", eta = 0.1),
      ret_extra = "P",
      epoch_callback = FALSE,
      verbose = FALSE
    )
    expect_equal(actual$Y, expected, tolerance = 1e-12)
    expect_identical(actual$P, initial$P)
    expect_identical(X, before$X)
    expect_identical(Y, before$Y)
  }
})

test_that("complete native MDS runs match independent distance updates", {
  n <- 31
  X <- matrix(sin(seq_len(n * 4)^2), n, 4)
  Y <- cbind(sin(seq_len(n)), cos(seq_len(n)))
  before <- unserialize(serialize(list(X = X, Y = Y), NULL))
  R <- unname(as.matrix(stats::dist(X)))
  expected <- Y
  eps <- 0.01
  for (step in 1:4) {
    D <- unname(as.matrix(stats::dist(expected)))
    G <- matrix(0, n, 2)
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        G[i, ] <- G[i, ] +
          4 *
            (D[i, j] - R[i, j]) /
            (D[i, j] + eps) *
            (expected[i, ] - expected[j, ])
      }
    }
    expected <- expected - 0.001 * G
    expected <- sweep(expected, 2, colMeans(expected))
  }
  for (threads in c(1L, 2L)) {
    actual <- smallvis(
      X,
      Y_init = Y,
      method = list("mmds", eps = eps),
      scale = FALSE,
      use_cpp = TRUE,
      n_threads = threads,
      max_iter = 4,
      opt = list("steepd", eta = 0.001),
      ret_extra = TRUE,
      epoch_callback = FALSE,
      verbose = FALSE
    )
    expect_equal(actual$Y, expected, tolerance = 1e-12)
    expect_equal(
      as.numeric(tail(actual$itercosts, 1)),
      sum((R - as.matrix(stats::dist(expected)))^2),
      tolerance = 1e-12
    )
    expect_identical(actual$iter, 4)
    expect_identical(X, before$X)
    expect_identical(Y, before$Y)
  }
})
