test_that("AB-SNE scalar objectives match independent non-self formulas", {
  X <- cbind(sin(1:8), cos(1:8))
  Y <- matrix(sin(1:16), 8)
  cases <- list(
    c(alpha = 0.75, lambda = 0.5),
    c(alpha = 1, lambda = 2),
    c(alpha = 1, lambda = 1),
    c(alpha = 2, lambda = 2),
    c(alpha = -0.5, lambda = -0.5),
    c(alpha = 1, lambda = 0),
    c(alpha = 2, lambda = 0),
    c(alpha = 0, lambda = 2),
    c(alpha = 0, lambda = -0.5),
    c(alpha = 0, lambda = 0)
  )
  for (parameters in cases) {
    result <- smallvis(
      X,
      method = c(list("absne"), as.list(parameters)),
      scale = FALSE,
      perplexity = 4,
      Y_init = Y,
      max_iter = 0,
      ret_extra = c("P", "Q"),
      epoch_callback = FALSE,
      verbose = FALSE
    )
    # All non-self entries stay above clipping thresholds, including P^2.
    off_diagonal <- row(result$P) != col(result$P)
    expect_gt(min(result$P[off_diagonal]^2), .Machine$double.eps)
    expect_gt(min(result$Q[off_diagonal]^2), .Machine$double.eps)
    expected <- absne_reference_cost(
      result$P,
      Y,
      parameters[["alpha"]],
      parameters[["lambda"]]
    )
    expect_equal(
      sum(result$costs),
      expected,
      tolerance = 1e-10,
      info = paste(parameters, collapse = ", ")
    )
    expect_gte(expected, 0)
  }
})

test_that("AB-SNE corrected branches have zero divergence when P equals Q", {
  for (parameters in list(
    c(alpha = 2, lambda = 2),
    c(alpha = -0.5, lambda = -0.5),
    c(alpha = 1, lambda = 0),
    c(alpha = 2, lambda = 0)
  )) {
    result <- smallvis(
      diag(8),
      method = c(list("absne"), as.list(parameters)),
      scale = FALSE,
      k = 8,
      perplexity = 7,
      Y_init = diag(8),
      max_iter = 0,
      ret_extra = c("P", "Q"),
      epoch_callback = FALSE,
      verbose = FALSE
    )
    expect_equal(result$P, result$Q, tolerance = 1e-14)
    expect_equal(result$costs, rep(0, 8), tolerance = 1e-10)
    expect_equal(
      absne_reference_cost(
        result$P,
        result$Y,
        parameters[["alpha"]],
        parameters[["lambda"]]
      ),
      0,
      tolerance = 1e-10
    )
  }
})
