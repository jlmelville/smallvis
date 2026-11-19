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

test_that("ABSSNE scores and gradients agree with the Gaussian AB divergence", {
  X <- cbind(sin(1:8), cos(1:8))
  Y <- matrix(sin(1:16), 8)
  direction <- Y / sqrt(sum(Y^2))
  cases <- rbind(
    c(1, 1),
    c(2, 2),
    c(0, 1),
    c(0, -0.5),
    c(1, 0),
    c(0, 0),
    c(0.75, 0.5)
  )
  for (i in seq_len(nrow(cases))) {
    alpha <- cases[i, 1]
    lambda <- cases[i, 2]
    result <- smallvis(
      X,
      method = list("abssne", alpha = alpha, lambda = lambda),
      scale = FALSE,
      perplexity = 4,
      Y_init = Y,
      max_iter = 0,
      ret_extra = "P",
      epoch_callback = FALSE,
      verbose = FALSE
    )
    ref <- function(coords) {
      absne_reference_cost(result$P, coords, alpha, lambda, gaussian = TRUE)
    }
    expect_equal(sum(result$costs), ref(Y), tolerance = 1e-10)
    cost <- cost_init(
      abssne(4, alpha = alpha, lambda = lambda),
      X,
      max_iter = 2
    )
    h <- 1e-5
    expected <- (ref(Y + h * direction) - ref(Y - h * direction)) / (2 * h)
    expect_equal(
      sum(cost_grad(cost, Y)$G * direction),
      expected,
      tolerance = 1e-7
    )
  }
})

test_that("Gaussian AB limits and the ABSNE threshold have zero self-divergence", {
  for (method in list(
    "abssne",
    list("abssne", alpha = 0, lambda = 0),
    list("absne", alpha = 1e-5, lambda = 1e-5)
  )) {
    result <- smallvis(
      diag(8),
      method = method,
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
  }
})
