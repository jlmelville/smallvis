test_that("coordinates and precomputed distances give the same Gaussian affinities", {
  X <- matrix(sin(seq_len(144)^2), 24, 6)
  run <- function(input, ...) {
    smallvis(
      input,
      scale = FALSE,
      perplexity = 4,
      max_iter = 0,
      Y_init = matrix(0, 24, 2),
      ret_extra = "P",
      epoch_callback = FALSE,
      verbose = FALSE,
      ...
    )
  }
  expect_equal(run(X)$P, run(stats::dist(X))$P, tolerance = 1e-10)
  expect_equal(
    run(stats::dist(X), use_cpp = TRUE)$P,
    run(stats::dist(X), use_cpp = FALSE)$P
  )
  expect_equal(
    run(X, method = list("tsne", inp_kernel = "exp"))$P,
    run(stats::dist(X), method = list("tsne", inp_kernel = "exp"))$P,
    tolerance = 1e-10
  )
})

test_that("distance kNN graphs select exactly k non-self neighbors in each row", {
  X <- matrix(sin(seq_len(144)^2), 24, 6)
  D <- as.matrix(stats::dist(X))
  diag(D) <- Inf
  for (k in c(1L, 4L)) {
    expected <- matrix(0, 24, 24)
    for (i in 1:24) {
      expected[i, order(D[i, ])[seq_len(k)]] <- 1
    }
    for (input in list(X, stats::dist(X))) {
      result <- smallvis(
        input,
        method = list("asne", inp_kernel = "knn"),
        scale = FALSE,
        perplexity = k,
        max_iter = 0,
        ret_extra = "P",
        epoch_callback = FALSE,
        verbose = FALSE
      )
      expect_equal(result$P, expected / k)
    }
    sparse <- knn_graph(stats::dist(X), k, ret_sparse = TRUE)
    expect_equal(as.matrix(sparse), expected)
  }
  ties <- stats::dist(matrix(0, 10, 2))
  graph <- knn_graph(ties, 4)
  expect_equal(diag(graph), rep(0, 10))
  expect_equal(rowSums(graph), rep(4, 10))
  expect_true(all(is.finite(graph)))
})

test_that("multiscale calibration preserves each method's affinity normalization", {
  X <- matrix(sin(seq_len(144)^2), 24, 6)
  for (method in c("tsne", "asne")) {
    run <- function(perplexity) {
      smallvis(
        X,
        method = method,
        scale = FALSE,
        perplexity = perplexity,
        Y_init = matrix(0, 24, 2),
        max_iter = 0,
        ret_extra = "P",
        epoch_callback = FALSE,
        verbose = FALSE
      )$P
    }
    expect_equal(run(4), run(list("multiscale", 4)), tolerance = 1e-12)
    expect_equal(
      run(list("multiscale", c(4, 6))),
      (run(4) + run(6)) / 2,
      tolerance = 1e-12
    )
    if (method == "asne") {
      expect_equal(rowSums(run(list("multiscale", 4))), rep(1, 24))
    } else {
      expect_equal(sum(run(list("multiscale", 4))), 1)
    }
  }
})

test_that("nearest-neighbor Gaussian calibration accepts distance inputs", {
  X <- matrix(sin(seq_len(144)^2), 24, 6)
  run <- function(input, bh = FALSE) {
    smallvis(
      input,
      method = list("tsne", inp_kernel = "nngaussian"),
      bh = bh,
      nn = "exact",
      perplexity = 4.5,
      scale = FALSE,
      max_iter = 0,
      ret_extra = "P",
      epoch_callback = FALSE,
      verbose = FALSE
    )$P
  }
  # rnndescent returns distances rounded to its native precision.
  coord_nn <- get_nn(X, k = 14, exact = TRUE)
  dist_nn <- get_nn(stats::dist(X), k = 14, exact = TRUE)
  expect_equal(coord_nn$idx, dist_nn$idx)
  expect_equal(coord_nn$dist, dist_nn$dist, tolerance = 1e-6)
  expect_equal(run(X), run(stats::dist(X)), tolerance = 1e-6)
  expect_equal(as.matrix(run(stats::dist(X), TRUE)), run(X), tolerance = 1e-6)
})
