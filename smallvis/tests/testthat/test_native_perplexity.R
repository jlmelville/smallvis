test_that("Gaussian perplexity is bounded by the non-self neighbor count", {
  X <- as.matrix(iris[1:12, 1:4])
  run <- function(perplexity = 30, use_cpp = TRUE) {
    smallvis(
      X,
      perplexity = perplexity,
      use_cpp = use_cpp,
      max_iter = 0,
      epoch_callback = FALSE,
      verbose = FALSE
    )
  }
  expect_error(run(), "Perplexity.*n - 1")
  expect_error(run(12), "Perplexity.*n - 1")
  expect_true(all(is.finite(run(11))))
  expect_error(run(c(rep(5, 11), 12), FALSE), "Perplexity.*n - 1")
  expect_error(run(list("multiscale", c(5, 12)), FALSE), "Perplexity.*n - 1")
  expect_error(run(list("idp", c(5, 12)), FALSE), "Perplexity.*n - 1")
  for (p in c(NA_real_, Inf, 0)) {
    expect_error(run(p), "Perplexity.*finite")
  }
  expect_true(all(is.finite(run(5))))

  # Registered calibrators must bound selection even without the R validator.
  expect_error(find_beta_cpp(X, 30), "Perplexity")
  expect_error(
    find_beta_knn_cpp(matrix(1, 12, 4), matrix(1L, 12, 4), 5),
    "Perplexity"
  )
})

test_that("native calibration failures retain finite non-self affinities", {
  X <- diag(24)
  dense <- find_beta_cpp(X, perplexity = 4)
  expect_equal(dense$bad_perp, 24L)
  expect_equal(diag(dense$W), rep(0, 24))
  expect_equal(rowSums(dense$W > 0), rep(4, 24))
  expect_true(all(is.finite(dense$W)))
  expect_equal(find_beta_cpp(X, 4, n_threads = 2), dense)

  # Tie-breaking must exclude self even when every input observation is equal.
  duplicates <- matrix(0, 24, 2)
  zero <- find_beta_cpp(duplicates, 4)
  expect_equal(zero$bad_perp, 24L)
  expect_equal(diag(zero$W), rep(0, 24))
  expect_equal(rowSums(zero$W > 0), rep(4, 24))
  expect_true(all(is.finite(zero$W)))

  indices <- t(vapply(
    seq_len(24),
    function(i) {
      as.integer(setdiff(seq_len(24), i)[seq_len(12)])
    },
    integer(12)
  ))
  for (distance in c(0, sqrt(2))) {
    knn <- find_beta_knn_cpp(matrix(distance, 24, 12), indices, 4)
    expect_equal(knn$bad_perp, 24L)
    expect_equal(
      knn,
      find_beta_knn_cpp(
        matrix(distance, 24, 12),
        indices,
        4,
        n_threads = 2
      )
    )
    expect_true(all(is.finite(knn$P)))
    expect_equal(diag(knn$P), rep(0, 24))
    expect_equal(rowSums(knn$P), rep(1, 24))
  }
})

test_that("calibration failures warn once and retain valid public affinities", {
  routes <- list(
    list("tsne", inp_kernel = "gaussian", use_cpp = FALSE),
    list("tsne", inp_kernel = "gaussian", use_cpp = TRUE),
    list("tsne", inp_kernel = "nngaussian", use_cpp = TRUE)
  )
  for (X in list(diag(24), matrix(0, 24, 4))) {
    for (perplexity in c(1, 4)) {
      for (method in routes) {
        warnings <- character()
        res <- withCallingHandlers(
          smallvis(
            X,
            Y_init = matrix(sin(seq_len(48)), 24),
            scale = FALSE,
            perplexity = perplexity,
            max_iter = 0,
            method = method,
            n_threads = 2,
            ret_extra = c("P", "V"),
            epoch_callback = FALSE,
            verbose = FALSE
          ),
          warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
          }
        )
        expect_equal(warnings, "24 perplexity calibrations failed")
        expect_equal(diag(res$P), rep(0, 24))
        expect_true(all(is.finite(res$P)))
        expect_true(all(rowSums(res$P) > 0))
        expect_equal(sum(res$P), 1)
        expect_true(all(is.finite(res$V)))
        expect_equal(diag(res$V), rep(0, 24))
        neighbors <- if (method$inp_kernel == "nngaussian") {
          3 * perplexity
        } else {
          perplexity
        }
        expect_equal(rowSums(res$V > 0), rep(neighbors, 24))
      }
    }
  }
})

test_that("dense Gaussian diagnostics retain calibrated precisions", {
  X <- matrix(sin(seq_len(96)), 24, 4)
  for (use_cpp in c(FALSE, TRUE)) {
    expect_warning(
      res <- smallvis(
        X,
        scale = FALSE,
        perplexity = 4,
        use_cpp = use_cpp,
        max_iter = 0,
        ret_extra = c("beta", "V"),
        epoch_callback = FALSE,
        verbose = FALSE
      ),
      NA
    )
    expect_length(res$beta, 24)
    expect_true(all(is.finite(res$beta)))
    D2 <- as.matrix(stats::dist(X))^2
    dimnames(D2) <- NULL
    W <- exp(-D2 * res$beta)
    diag(W) <- 0
    expect_equal(res$V, W, tolerance = 1e-12)
    expect_equal(perpp(W), rep(4, 24), tolerance = 1e-4)
  }
})

test_that("native kNN output checks the matrix shape and write indices", {
  # These checks protect the native scatter into an n-by-n affinity matrix.
  expect_error(
    find_beta_knn_cpp(matrix(1, 12, 4), matrix(1L, 12, 3), 3),
    "same dimensions"
  )
  for (index in c(0L, NA_integer_, 13L)) {
    expect_error(
      find_beta_knn_cpp(matrix(1, 12, 4), matrix(index, 12, 4), 3),
      "out of range"
    )
  }
})

test_that("fixed bandwidths retain their own parameter domain", {
  res <- smallvis(
    iris10,
    method = list("tsne", inp_kernel = "sigma"),
    perplexity = 100,
    use_cpp = TRUE,
    max_iter = 0,
    ret_extra = "P",
    epoch_callback = FALSE,
    verbose = FALSE
  )
  expect_true(all(is.finite(res$P)))
  expect_equal(sum(res$P), 1)
})
