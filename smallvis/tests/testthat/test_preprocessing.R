test_that("scaling treats matrices and data frames as columns", {
  X <- cbind(a = 1:24, b = cos(1:24), constant = 7)
  expected <- scale(X[, 1:2])
  dimnames(expected) <- NULL
  for (input in list(X, as.data.frame(X))) {
    res <- smallvis(
      input,
      scale = TRUE,
      perplexity = 4,
      max_iter = 0,
      ret_extra = "X",
      epoch_callback = FALSE,
      verbose = FALSE
    )
    actual <- res$X
    dimnames(actual) <- NULL
    attr(actual, "scaled:center") <- attr(actual, "scaled:scale") <- NULL
    attr(expected, "scaled:center") <- attr(expected, "scaled:scale") <- NULL
    expect_equal(actual, expected)
  }
  one <- smallvis(
    data.frame(a = 1:24, label = letters[1:24]),
    k = 1,
    scale = TRUE,
    perplexity = 4,
    max_iter = 0,
    ret_extra = "X",
    epoch_callback = FALSE,
    verbose = FALSE
  )
  expect_equal(dim(one$X), c(24L, 1L))
})

test_that("absmax scaling uses the largest absolute centered value", {
  X <- matrix(c(-9, rep(1, 9)), 10, 1)
  res <- smallvis(
    X,
    scale = "absmax",
    k = 1,
    perplexity = 9,
    max_iter = 0,
    ret_extra = "X",
    epoch_callback = FALSE,
    verbose = FALSE
  )
  expect_equal(as.numeric(res$X), as.numeric(X / 9))
  for (scale in list(TRUE, "absmax", "range")) {
    expect_error(
      smallvis(
        matrix(7, 24, 2),
        scale = scale,
        perplexity = 4,
        max_iter = 0,
        epoch_callback = FALSE,
        verbose = FALSE
      ),
      "Scaling|scaling"
    )
  }
})

test_that("PCA metadata and initialization reflect whether preprocessing ran", {
  X <- matrix(sin(seq_len(96)^2), 24, 4)
  res <- smallvis(
    X,
    scale = FALSE,
    pca = TRUE,
    initial_dims = 50,
    Y_init = "pca",
    perplexity = 4,
    max_iter = 0,
    ret_extra = TRUE,
    epoch_callback = FALSE,
    verbose = FALSE
  )
  expect_null(res$pca_dims)
  expect_equal(
    as.matrix(stats::dist(res$Y)),
    as.matrix(stats::dist(stats::prcomp(X)$x[, 1:2])),
    ignore_attr = TRUE,
    tolerance = 1e-10
  )
  for (k in c(1L, 3L)) {
    res <- smallvis(
      X,
      scale = FALSE,
      pca = TRUE,
      initial_dims = 3,
      Y_init = "pca",
      k = k,
      perplexity = 4,
      max_iter = 0,
      ret_extra = c("X"),
      epoch_callback = FALSE,
      verbose = FALSE
    )
    expect_equal(dim(res$Y), c(24L, k))
    expect_equal(res$pca_dims, 3)
    expect_equal(res$origD, 4)
    expect_equal(
      res$Y,
      unname(res$X[, seq_len(k), drop = FALSE]),
      ignore_attr = TRUE
    )
  }
})

test_that("truncated PCA whitening matches full SVD eigenvalues and scores", {
  X <- matrix(sin(seq_len(24 * 14)^2), 24, 14)
  set.seed(42)
  res <- smallvis(
    X,
    scale = FALSE,
    pca = "whiten",
    initial_dims = 3,
    Y_init = "pca",
    k = 3,
    perplexity = 4,
    max_iter = 0,
    ret_extra = "X",
    epoch_callback = FALSE,
    verbose = FALSE
  )
  full <- stats::prcomp(X)
  lambda <- full$sdev[1:3]^2
  expected <- sweep(full$x[, 1:3], 2, sqrt(lambda + 1e-5), "/")
  expect_equal(res$whiten_dims, 3)
  expect_equal(res$origD, 14)
  expect_equal(
    as.matrix(stats::dist(res$X)),
    as.matrix(stats::dist(expected)),
    ignore_attr = TRUE,
    tolerance = 1e-7
  )
  expect_equal(
    unname(diag(stats::cov(res$X))),
    lambda / (lambda + 1e-5),
    tolerance = 1e-7
  )
  one <- pca_scores(X[, 1:2], ncol = 1)
  expect_equal(dim(one), c(24L, 1L))
})

test_that("PCA initialization rejects unavailable components", {
  expect_error(
    smallvis(
      matrix(1:24, 24, 1),
      Y_init = "pca",
      k = 2,
      perplexity = 4,
      max_iter = 0,
      epoch_callback = FALSE,
      verbose = FALSE
    ),
    "between 1 and 1 components"
  )
})
