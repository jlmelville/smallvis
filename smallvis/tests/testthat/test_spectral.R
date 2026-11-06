test_that("normalized laplacian", {
  # These numbers come from running UMAP Python code:
  # spectral_layout(pairwise_distances(iris.data[0:10, :]))
  # NB:
  # 1. iris data in scikit-learn is currently from UCI repo, which has errors
  #   (although this doesn't affect the first ten entries)
  # 2. eigenvector calculation is not that converged and specifies a starting
  #   vector that we can't supply with either RSpectra or eigen.
  # 3. The eigenvectors are only identical up to a sign, so we take the absolute
  #   values.
  # fmt: skip
  expected_norm_lap <-
  c2y(0.7477, -0.1292, -0.03001, 0.02127, -0.563, -0.01149, 0.1402,
    -0.2725, -0.01241, 0.1084, -0.106, -0.5723, 0.2024, -0.3082,
    0.1642, -5.549e-05, -0.04843, -0.1747, 0.1684, 0.6611)

  res <- normalized_spectral_init(x2d(iris[1:10, ]), use_RSpectra = TRUE)
  expect_equal(abs(res), abs(expected_norm_lap), tolerance = 1e-2)

  res <- normalized_spectral_init(x2d(iris[1:10, ]), use_RSpectra = FALSE)
  expect_equal(abs(res), abs(expected_norm_lap), tolerance = 1e-2)
})

test_that("laplacian eigenmap uses the leading algebraic modes", {
  A <- x2d(iris[1:10, ])
  degrees <- colSums(A)
  values <- sort(Re(eigen(A / degrees)$values), decreasing = TRUE)[2:3]
  for (use_RSpectra in c(TRUE, FALSE)) {
    Y <- laplacian_eigenmap(A, use_RSpectra = use_RSpectra)
    expect_equal(A %*% Y, sweep(Y * degrees, 2, values, "*"), tolerance = 1e-9)
    expect_equal(colSums(Y * degrees), c(0, 0), tolerance = 1e-9)
  }
})

test_that("spectral initialization accepts sparse affinity matrices", {
  affinity <- Matrix::Matrix(x2d(iris[1:10, ]), sparse = TRUE)

  for (init in list(laplacian_eigenmap, normalized_spectral_init)) {
    dense <- init(as.matrix(affinity), use_RSpectra = FALSE)
    sparse <- init(affinity, use_RSpectra = FALSE)

    expect_equal(abs(sparse), abs(dense), tolerance = 1e-10)
  }
})
