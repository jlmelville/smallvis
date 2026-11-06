test_that("Laplacian eigenmaps select the smooth cycle eigenspace", {
  n <- 10
  A <- matrix(0, n, n)
  A[cbind(1:n, c(2:n, 1))] <- 1
  A <- A + t(A)
  angle <- 2 * pi * (0:9) / n
  basis <- qr.Q(qr(cbind(sin(angle), cos(angle))))
  projector <- tcrossprod(basis)
  for (use_RSpectra in c(FALSE, TRUE)) {
    Y <- laplacian_eigenmap(A, use_RSpectra = use_RSpectra)
    expect_equal(tcrossprod(qr.Q(qr(Y))), projector, tolerance = 1e-10)
    expect_equal(A %*% Y / 2, Y * cos(2 * pi / n), tolerance = 1e-10)
    expect_equal(colSums(Y), c(0, 0), tolerance = 1e-10)
  }
})

test_that("Laplacian eigenmaps satisfy the generalized problem with unequal degrees", {
  A <- outer(1:10, 1:10, function(i, j) exp(-abs(i - j) / 3))
  diag(A) <- 0
  degrees <- colSums(A)
  expect_gt(diff(range(degrees)), 0.1)
  eigenvalues <- sort(Re(eigen(A / degrees)$values), decreasing = TRUE)[2:3]
  for (use_RSpectra in c(FALSE, TRUE)) {
    Y <- laplacian_eigenmap(A, use_RSpectra = use_RSpectra)
    expect_equal(
      (diag(degrees) - A) %*% Y,
      sweep(Y * degrees, 2, 1 - eigenvalues, "*"),
      tolerance = 1e-10
    )
    expect_equal(colSums(Y * degrees), c(0, 0), tolerance = 1e-10)
    expect_equal(unname(crossprod(Y, Y * degrees)[1, 2]), 0, tolerance = 1e-10)
  }
})

test_that("spectral layouts are independent of requested optional diagnostics", {
  X <- matrix(sin(seq_len(144)^2), 24, 6)
  for (method in c("tsne", "asne")) {
    for (init in c("laplacian", "normlaplacian")) {
      run <- function(extra) {
        set.seed(42)
        smallvis(
          X,
          method = method,
          perplexity = 4,
          Y_init = init,
          max_iter = 0,
          ret_extra = extra,
          epoch_callback = FALSE,
          verbose = FALSE
        )$Y
      }
      expect_equal(
        as.matrix(stats::dist(run(TRUE))),
        as.matrix(stats::dist(run("V"))),
        tolerance = 1e-10
      )
    }
  }
})

test_that("small spectral inputs use a valid dense eigensolve", {
  A <- matrix(1, 3, 3) - diag(3)
  expect_equal(dim(laplacian_eigenmap(A, ndim = 2)), c(3L, 2L))
})

test_that("small normalized spectral inputs use a valid dense eigensolve", {
  A <- matrix(1, 3, 3) - diag(3)
  expect_equal(dim(normalized_spectral_init(A, ndim = 2)), c(3L, 2L))
})
