test_that("squared self-distances are exactly zero", {
  # Matrix multiplication can leave a small positive diagonal on some BLAS
  # implementations. Self-distances are structurally zero, not approximate.
  Y <- iris10_Y
  Y[2, 1] <- Y[2, 1] + .Machine$double.eps ^ (1 / 3)

  expect_equal(diag(safe_dist2(Y)), rep(0, nrow(Y)), tolerance = 0)
})
