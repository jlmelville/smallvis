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

test_that("EE and g-SNE construct equivalent targets from coordinates and distances", {
  set.seed(42)
  X <- matrix(rnorm(24 * 4), 24)
  Y <- matrix(rnorm(48), 24)
  D <- as.matrix(stats::dist(X))
  dimnames(D) <- NULL
  targets <- list(ee = D / sum(D), gsne = 1 + D^2)
  diag(targets$gsne) <- 0
  targets$gsne <- targets$gsne / sum(targets$gsne)
  for (method in c("ee", "gsne")) {
    field <- if (method == "ee") "Vn" else "Phat"
    run <- function(input, iterations = 0) {
      gradient <- NULL
      target <- NULL
      expect_warning(
        fit <- smallvis(
          input,
          method = method,
          scale = FALSE,
          perplexity = 4,
          Y_init = Y,
          use_cpp = FALSE,
          n_threads = 1,
          max_iter = iterations,
          opt = list("adam", eta = 0.001),
          epoch = 1,
          tol = 0,
          min_cost = -Inf,
          ret_extra = "P",
          epoch_callback = function(Y, iter, cost, cost_fn) {
            if (iter == 0) {
              target <<- cost_fn[[field]]
              gradient <<- cost_fn$gr(cost_fn, Y)$G
            }
          },
          verbose = FALSE
        ),
        NA
      )
      list(fit = fit, gradient = gradient, target = target)
    }
    coordinates <- run(X)
    distances <- run(stats::dist(X))
    for (result in list(coordinates, distances)) {
      expect_equal(dim(result$target), c(24L, 24L))
      expect_equal(
        result$target,
        targets[[method]],
        ignore_attr = TRUE,
        tolerance = 1e-12
      )
    }
    expect_equal(coordinates$fit$P, distances$fit$P, tolerance = 1e-10)
    expect_equal(
      coordinates$fit$costs,
      distances$fit$costs,
      ignore_attr = TRUE,
      tolerance = 1e-10
    )
    expect_equal(
      coordinates$gradient,
      distances$gradient,
      ignore_attr = TRUE,
      tolerance = 1e-10
    )
    coordinate_run <- run(X, 4)$fit
    distance_run <- run(stats::dist(X), 4)$fit
    expect_equal(coordinate_run$iter, 4)
    expect_equal(distance_run$iter, 4)
    expect_equal(
      coordinate_run$Y,
      distance_run$Y,
      ignore_attr = TRUE,
      tolerance = 1e-10
    )
    expect_equal(
      coordinate_run$itercosts,
      distance_run$itercosts,
      tolerance = 1e-10
    )
  }
})

test_that("geodesic MDS supports one non-self neighbor with either input", {
  X <- cbind((1:8)^2, 0)
  expected <- as.matrix(stats::dist(X))
  # Increasing gaps make the undirected one-neighbor graph a connected chain.
  for (input in list(X, stats::dist(X))) {
    fit <- smallvis(
      input,
      method = list("gmmds", k = 1),
      scale = FALSE,
      Y_init = cbind(sin(1:8), cos(1:8)),
      max_iter = 2,
      opt = list("adam", eta = 0.001),
      epoch = 1,
      tol = 0,
      ret_extra = "geo",
      epoch_callback = FALSE,
      verbose = FALSE
    )
    expect_equal(fit$geo, expected, ignore_attr = TRUE, tolerance = 1e-12)
    expect_equal(fit$iter, 2)
    expect_true(all(is.finite(fit$Y)))
    expect_true(all(is.finite(fit$costs)))
  }
})
