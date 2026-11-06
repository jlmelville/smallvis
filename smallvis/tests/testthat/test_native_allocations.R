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
