library(smallvis)
context("API")

expected_Y <- matrix(c(
  -188.6, 160.2, 81.5, 182.8, -247, -266.6, 41.8, -125.2, 298.6, 62.5,
  -25.6, 62.1, -66.6, -94.2, 54, 190.9, -166.5, 61.4, -93.1, 77.7
), ncol = 2)

test_that("basic return value", {
i10_tsne <- smallvis(iris10, Y_init = iris10_Y, perplexity = 5,
                     epoch_callback = NULL, verbose = FALSE)
  expect_equal(i10_tsne, expected_Y, tolerance = 0.1)
})


test_that("extra return values", {
  i10_tsne <- smallvis(iris10, Y_init = iris10_Y, perplexity = 5,
                       epoch_callback = NULL, verbose = FALSE,
                       ret_extra = TRUE)
  expect_equal(i10_tsne$N, 10)
  expect_equal(i10_tsne$origD, 4)
  expect_equal(i10_tsne$scale, "absmax")
  expect_equal(i10_tsne$Y_init, "matrix")
  expect_equal(i10_tsne$method, "tsne")
  expect_equal(i10_tsne$perplexity, 5)
  expect_equal(i10_tsne$stop_lying_iter, 100)
  expect_equal(i10_tsne$exaggeration_factor, 1)

  expected_itercosts <- c(0.3096, 0.05035, 0.04133, 0.04129, 0.04129)
  names(expected_itercosts) <- c(100, 200, 300, 400, 500)
  expect_equal(i10_tsne$itercosts, expected_itercosts, tolerance = 1e-5)

  expect_equal(i10_tsne$costs, c(0.005045, 0.0005255, -0.005447, 0.003863,
                                 0.0005298,  0.006410,  0.01191, 0.007799,
                                 0.003098,  0.007556), tolerance = 1e-5)

  expect_equal(i10_tsne$Y, expected_Y, tolerance = 0.1)
})
