# The fixed symmetric P isolates tree calculations from neighbor calibration.
bh_probabilities <- function(n) {
  P <- outer(seq_len(n), seq_len(n), "+")
  diag(P) <- 0
  P / sum(P)
}

bh_dense_reference <- function(P, Y, eps = .Machine$double.eps) {
  W <- 1 / (1 + as.matrix(stats::dist(Y))^2)
  dimnames(W) <- NULL
  diag(W) <- 0
  Q <- W / (sum(W) + eps)
  used <- P > 0
  G <- repulsion <- matrix(0, nrow(Y), 2)
  for (i in seq_len(nrow(Y))) {
    direction <- matrix(Y[i, ], nrow(Y), 2, byrow = TRUE) - Y
    repulsion[i, ] <- colSums(direction * W[i, ]^2) / (sum(W) + eps)
    G[i, ] <- 4 * colSums(direction * (P[i, ] - Q[i, ]) * W[i, ])
  }
  list(
    cost = sum(P[used] * log(P[used] / Q[used])),
    gradient = G,
    repulsion = repulsion
  )
}

bh_native_result <- function(
  P,
  Y,
  eps = .Machine$double.eps,
  n_threads = 1,
  theta = 0
) {
  sparse <- as(Matrix::Matrix(P, sparse = TRUE), "generalMatrix")
  args <- list(
    sparse@i,
    sparse@p,
    sparse@x,
    Y,
    theta = theta,
    eps = eps,
    n_threads = n_threads
  )
  list(
    cost = sum(P[P > 0] * log(P[P > 0])) - do.call(bh_plogq_cpp, args),
    gradient = 4 * do.call(bh_tsne_gradient_cpp, args)
  )
}

expect_bh_reference <- function(Y, orders = list(seq_len(nrow(Y)))) {
  P <- bh_probabilities(nrow(Y))
  expected <- bh_dense_reference(P, Y)
  for (order in orders) {
    for (threads in c(1L, 2L)) {
      actual <- bh_native_result(
        P[order, order],
        Y[order, ],
        n_threads = threads
      )
      expect_equal(actual$cost, expected$cost, tolerance = 1e-12)
      expect_equal(
        actual$gradient[order(order), ],
        expected$gradient,
        tolerance = 1e-12
      )
    }
  }
}

test_that("coincident leaves exclude exactly one query observation", {
  Y <- matrix(0, 24, 2)
  expect_bh_reference(Y)
  res <- smallvis(
    as.matrix(iris[1:24, 1:4]),
    Y_init = Y,
    perplexity = 4,
    bh = TRUE,
    nn = "exact",
    theta = 0,
    max_iter = 0,
    ret_extra = "P",
    epoch_callback = FALSE,
    verbose = FALSE
  )
  P <- as.matrix(res$P)
  used <- P > 0
  expected <- sum(P[used] * log(P[used] * 24 * 23))
  expect_equal(res$costs, expected, tolerance = 1e-12)
  expect_equal(bh_native_result(P, Y)$gradient, Y, tolerance = 1e-14)
})

test_that("summary nodes containing the query are expanded", {
  Y <- rbind(c(0, 0), c(10, 0))
  P <- bh_probabilities(2)
  actual <- bh_native_result(P, Y, theta = 10)
  expected <- bh_dense_reference(P, Y)
  expect_equal(actual$cost, expected$cost, tolerance = 1e-12)
  expect_equal(actual$gradient, expected$gradient, tolerance = 1e-12)
})

test_that("splitting a leaf preserves all duplicate observations", {
  Y <- rbind(c(0, 0), c(0, 0), c(0, 0), c(2, 1), c(-3, 4), c(1, -2))
  expect_bh_reference(Y, list(1:6, 6:1, c(4, 1, 6, 2, 5, 3)))
})

test_that("points at internal centroids still enter the children", {
  Y <- rbind(c(-1, 0), c(1, 0), c(0, 0), c(3, 1))
  expect_bh_reference(Y, list(1:4, c(4, 3, 1, 2)))
  control <- cbind(sin(1:12), cos(1:12))
  expect_bh_reference(control, list(1:12, 12:1))
})

test_that("Barnes-Hut epsilon stabilizes the normalizer without merging points", {
  Y <- cbind(sin(1:24), cos(1:24)) * 1e-4
  P <- bh_probabilities(24)
  for (eps in c(.Machine$double.eps, 1e-3)) {
    actual <- bh_native_result(P, Y, eps)
    expected <- bh_dense_reference(P, Y, eps)
    # Relative comparison after removing attraction detects a missing tiny force.
    W <- 1 / (1 + as.matrix(stats::dist(Y))^2)
    dimnames(W) <- NULL
    diag(W) <- 0
    K <- P * W
    attraction <- Y * rowSums(K) - K %*% Y
    repulsion <- attraction - actual$gradient / 4
    expect_equal(
      repulsion / max(abs(expected$repulsion)),
      expected$repulsion / max(abs(expected$repulsion)),
      tolerance = 1e-12
    )
    expect_equal(actual$cost, expected$cost, tolerance = 1e-12)
    expect_equal(actual$gradient, expected$gradient, tolerance = 1e-12)

    res <- smallvis(
      as.matrix(iris[1:24, 1:4]),
      Y_init = Y,
      bh = TRUE,
      nn = "exact",
      perplexity = 4,
      theta = 0,
      eps = eps,
      max_iter = 0,
      ret_extra = "P",
      epoch_callback = FALSE,
      verbose = FALSE
    )
    expect_equal(
      res$costs,
      bh_dense_reference(as.matrix(res$P), Y, eps)$cost,
      tolerance = 1e-12
    )
    finite_difference <- function(step) {
      forward <- backward <- Y
      forward[1, 1] <- forward[1, 1] + step
      backward[1, 1] <- backward[1, 1] - step
      (bh_native_result(P, forward, eps)$cost -
        bh_native_result(P, backward, eps)$cost) /
        (2 * step)
    }
    for (step in c(1e-4, 1e-5)) {
      expect_equal_abs(
        finite_difference(step),
        actual$gradient[1, 1],
        tolerance = 1e-8
      )
    }
  }
})

test_that("Barnes-Hut rejects non-finite coordinates in a bounded subprocess", {
  dev <- pkgload::is_dev_package("smallvis")
  path <- system.file(package = "smallvis")
  result <- callr::r(
    function(dev, path) {
      if (dev) {
        pkgload::load_all(path, quiet = TRUE)
      } else {
        library("smallvis", lib.loc = dirname(path))
      }
      stopifnot(
        normalizePath(system.file(package = "smallvis")) == normalizePath(path)
      )
      P <- outer(1:4, 1:4, "+")
      diag(P) <- 0
      sparse <- as(Matrix::Matrix(P / sum(P), sparse = TRUE), "generalMatrix")
      errors <- list(
        public = character(),
        cost = character(),
        gradient = character()
      )
      capture_error <- function(expr) {
        tryCatch(
          {
            force(expr)
            "returned invalid coordinates"
          },
          error = conditionMessage
        )
      }
      for (bad in c(NA_real_, NaN, Inf, -Inf)) {
        for (column in 1:2) {
          Y <- cbind(sin(1:24), cos(1:24))
          Y[2, column] <- bad
          errors$public <- c(
            errors$public,
            capture_error(smallvis::smallvis(
              as.matrix(iris[1:24, 1:4]),
              Y_init = Y,
              bh = TRUE,
              nn = "exact",
              perplexity = 4,
              max_iter = 0,
              epoch_callback = FALSE,
              verbose = FALSE
            ))
          )
          args <- list(sparse@i, sparse@p, sparse@x, Y[1:4, ])
          errors$cost <- c(
            errors$cost,
            capture_error(do.call(
              smallvis:::bh_plogq_cpp,
              args
            ))
          )
          errors$gradient <- c(
            errors$gradient,
            capture_error(do.call(
              smallvis:::bh_tsne_gradient_cpp,
              args
            ))
          )
        }
      }
      errors
    },
    args = list(dev = dev, path = path),
    libpath = .libPaths(),
    timeout = 30
  )
  for (route in c("public", "cost", "gradient")) {
    expect_equal(
      result[[route]],
      rep("Barnes-Hut coordinates must be finite", 8)
    )
  }
})

test_that("Barnes-Hut defaults to theta 0.5 and returns finite results", {
  run <- function(...) {
    set.seed(42)
    smallvis(
      iris,
      bh = TRUE,
      nn = "exact",
      perplexity = 10,
      Y_init = "spca",
      max_iter = 10,
      n_threads = 2,
      ret_extra = TRUE,
      epoch_callback = FALSE,
      verbose = FALSE,
      ...
    )
  }
  result <- run()
  expect_identical(formals(smallvis)$theta, 0.5)
  expect_identical(result$Y, run(theta = 0.5)$Y)
  expect_true(all(is.finite(result$Y)))
  expect_true(all(is.finite(result$itercosts)))
})

test_that("finite-resolution leaves finish and retain exact pair contributions", {
  dev <- pkgload::is_dev_package("smallvis")
  path <- system.file(package = "smallvis")
  # Keep the recursion regression bounded even if tree construction breaks again.
  results <- callr::r(
    function(dev, path) {
      if (dev) {
        pkgload::load_all(path, quiet = TRUE)
      } else {
        library("smallvis", lib.loc = dirname(path))
      }
      eps <- .Machine$double.eps
      Y <- cbind(c(rep(1 + eps, 23), 1 + 2 * eps), 0)
      P <- outer(1:24, 1:24, "+")
      diag(P) <- 0
      P <- P / sum(P)
      results <- list()
      for (order in list(1:24, 24:1, c(24, 1:23))) {
        sparse <- as(
          Matrix::Matrix(P[order, order], sparse = TRUE),
          "generalMatrix"
        )
        for (threads in c(1L, 2L)) {
          args <- list(
            sparse@i,
            sparse@p,
            sparse@x,
            Y[order, ],
            theta = 0,
            eps = eps,
            n_threads = threads
          )
          gradient <- 4 * do.call(smallvis:::bh_tsne_gradient_cpp, args)
          results[[length(results) + 1L]] <- list(
            cost = sum(P[P > 0] * log(P[P > 0])) -
              do.call(smallvis:::bh_plogq_cpp, args),
            gradient = gradient[order(order), ]
          )
        }
      }
      public <- smallvis::smallvis(
        cbind(1:24, sin(1:24)),
        scale = FALSE,
        Y_init = Y,
        perplexity = 4,
        bh = TRUE,
        nn = "exact",
        theta = 0,
        max_iter = 0,
        ret_extra = "P",
        epoch_callback = FALSE,
        verbose = FALSE
      )
      list(Y = Y, P = P, native = results, public = public)
    },
    args = list(dev = dev, path = path),
    libpath = .libPaths(),
    timeout = 30
  )
  expected <- bh_dense_reference(results$P, results$Y)
  for (actual in results$native) {
    expect_equal(actual$cost, expected$cost, tolerance = 1e-12)
    # The forces are of order machine epsilon: an absolute bound would admit zero.
    force_scale <- max(abs(expected$gradient))
    expect_gt(force_scale, 0)
    expect_equal(
      actual$gradient / force_scale,
      expected$gradient / force_scale,
      tolerance = 1e-12
    )
  }
  expect_equal(
    results$public$costs,
    bh_dense_reference(as.matrix(results$public$P), results$Y)$cost,
    tolerance = 1e-12
  )
})

test_that("new Barnes-Hut methods match dense updates on the same sparse graph", {
  X <- cbind(sin(1:24), cos(1:24), seq_len(24) / 7)
  Y <- cbind(sin(1:24) / 2, cos(1:24) / 3)
  Y[2, ] <- Y[1, ]
  for (method in list(list("tee", lambda = 0.04), list("hssne", alpha = 0.7))) {
    for (kernel in c("knn", "nngaussian")) {
      fit <- function(bh, theta = 0) {
        smallvis(
          X,
          method = method,
          bh = bh,
          inp_kernel = kernel,
          nn = "exact",
          perplexity = 4,
          Y_init = Y,
          max_iter = 1,
          eta = 0.5,
          theta = theta,
          # Expose an incorrect input floor in the reported objective.
          eps = 1e-5,
          n_threads = 2,
          ret_extra = "P",
          epoch_callback = FALSE,
          verbose = FALSE
        )
      }
      dense <- fit(FALSE)
      exact <- fit(TRUE)
      expect_s4_class(exact$P, "dgCMatrix")
      expect_equal(as.matrix(exact$P), dense$P, tolerance = 1e-12)
      expect_equal(exact$Y, dense$Y, tolerance = 1e-12)
      expect_equal(exact$costs, sum(dense$costs), tolerance = 1e-12)
      if (method[[1L]] == "tee") {
        P <- as.matrix(exact$P)
        W <- 1 / (1 + as.matrix(stats::dist(exact$Y))^2)
        diag(W) <- 0
        Q <- method$lambda * W
        used <- P > 0
        expected <- (sum(P[used] * log(P[used] / Q[used])) +
          sum(Q) -
          sum(P)) /
          sum(P)
        expect_equal(exact$costs, expected, tolerance = 1e-12)
      }
      approximate <- fit(TRUE, theta = 0.7)
      expect_true(all(is.finite(approximate$Y)))
      expect_true(is.finite(approximate$costs))
    }
  }
})

test_that("TEE has zero cost when scaled output weights equal input weights", {
  X <- cbind(0:1, 0)
  for (bh in c(FALSE, TRUE)) {
    run <- function(lambda) {
      smallvis(
        X,
        method = list("tee", lambda = lambda),
        bh = bh,
        inp_kernel = "knn",
        perplexity = 1,
        Y_init = X,
        theta = 0,
        max_iter = 0,
        ret_extra = TRUE,
        epoch_callback = FALSE,
        verbose = FALSE
      )
    }
    expect_equal(sum(run(2)$costs), 0, tolerance = 1e-14)
    expect_error(run(0), "TEE requires a finite positive scalar lambda")
    expect_error(run(-1), "TEE requires a finite positive scalar lambda")
  }
})

test_that("new Barnes-Hut methods require a sparse input kernel and two dimensions", {
  for (method in c("tee", "hssne")) {
    expect_error(
      smallvis(iris, method = method, bh = TRUE, inp_kernel = "gaussian"),
      "nngaussian"
    )
    expect_error(
      smallvis(iris, method = method, bh = TRUE, k = 3),
      "Barnes-Hut requires k = 2"
    )
  }
})
