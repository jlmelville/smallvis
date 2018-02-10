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

test_that("mmds", {
  res <- smallvis(iris10, Y_init = iris10_Y, method = "mmds", eta = 0.1,
                       epoch_callback = NULL, verbose = FALSE, ret_extra = TRUE)
  expect_equal(res$Y, c2y(0.5031, -0.3515, -0.3407, -0.59, 0.4913, 1.463, -0.209, 0.2736,
                            -1.056, -0.1838, -0.02805, -0.4261, 0.1157, 0.01604, 0.2117,
                            0.0396, 0.4632, -0.09042, 0.0585, -0.3602), tolerance = 1e-4)
  expect_equal(final_cost(res), 0.1715, tolerance = 1e-5)
})

test_that("sammon", {
  res <- smallvis(iris10, Y_init = iris10_Y, method = "sammon", eta = 0.1,
                  epoch_callback = NULL, verbose = FALSE, ret_extra = TRUE)
  expect_equal(res$Y, c2y(-0.5143, 0.3886, 0.3202, 0.6053, -0.4982, -1.473, 0.2011, -0.2667,
                            1.061, 0.1761, 0.02941, 0.4171, -0.1152, -0.0331, -0.2058, 0.03922,
                            -0.4762, 0.09549, -0.1027, 0.3518), tolerance = 1e-3)
  expect_equal(final_cost(res), 0.003427, tolerance = 1e-5)
})

test_that("gmmds", {
  res <- smallvis(iris10, Y_init = iris10_Y, method = "gmmds", eta = 0.1,
                  perplexity = 3,
                  epoch_callback = NULL, verbose = FALSE, ret_extra = TRUE)
  expect_equal(res$Y, c2y(0.07575, 0.3989, -0.02977, -0.1897, -0.2063, 0.09668, -0.484,
                            0.07552, -0.0509, 0.3138, 0.6367, -0.4354, -0.5816, -0.6774,
                            0.6217, 1.604, -0.1778, 0.386, -1.201, -0.1749), tolerance = 1e-4)
  expect_equal(final_cost(res), 0.4881, tolerance = 1e-4)
})

test_that("asne", {
  res <- smallvis(iris10, Y_init = iris10_Y, method = "asne", eta = 0.1,
                  perplexity = 5,
                  epoch_callback = NULL, verbose = FALSE, ret_extra = TRUE)
  expect_equal(res$Y, c2y(-0.9093, 0.5784, 0.4561, 0.7983, -0.9089, -1.669, 0.1867, -0.615,
                            1.695, 0.3876, 0.3947, 0.6318, -0.2393, -0.2071, 0.3158, -0.9542,
                            -0.683, 0.3709, -0.2324, 0.6029), tolerance = 1e-4)
  expect_equal(final_cost(res), 0.5177, tolerance = 1e-4)
})

test_that("ssne", {
  res <- smallvis(iris10, Y_init = iris10_Y, method = "ssne", eta = 0.1,
                  perplexity = 5,
                  epoch_callback = NULL, verbose = FALSE, ret_extra = TRUE)
  expect_equal(res$Y, c2y(-0.9462, 0.536, 0.4726, 0.8666, -0.9534, -1.266, 0.1875, -0.6616,
                        1.42, 0.3446, 0.2538, 0.6339, -0.2935, -0.2042, 0.06796, -0.5964,
                        -0.6855, 0.3511, -0.1556, 0.6284), tolerance = 1e-4)
  expect_equal(final_cost(res), 0.01683, tolerance = 1e-5)
})

test_that("hssne", {
  res <- smallvis(iris10, Y_init = iris10_Y, method = "hssne", eta = 0.1,
                  perplexity = 5,
                  epoch_callback = NULL, verbose = FALSE, ret_extra = TRUE)
  expect_equal(res$Y, c2y(-1.439, 0.8917, 0.793, 1.404, -1.35, -2.493, 0.4382, -0.9064,
                            2.204, 0.457, 0.2329, 0.8589, -0.3616, -0.226, -0.3591, -0.1309,
                            -0.9737, 0.3179, -0.1616, 0.8032), tolerance = 1e-3)
  expect_equal(final_cost(res), 0.01396, tolerance = 1e-5)
})

test_that("dhssne", {
  res <- smallvis(iris10, Y_init = iris10_Y, method = "dhssne", eta = 0.1,
                  perplexity = 5, epoch = 50,
                  epoch_callback = NULL, verbose = FALSE, ret_extra = c("alpha"))
  expect_equal(res$Y, c2y(-1.278, 0.784, 0.725, 1.279, -1.224, -2.276, 0.3911, -0.8396,
                        2.016, 0.4235, 0.1702, 0.7994, -0.3464, -0.2099, -0.3509, -0.07586,
                        -0.9089, 0.3002, -0.1318, 0.7539), tolerance = 1e-3)
  expect_equal(final_cost(res), 0.01270, tolerance = 1e-5)
  expect_equal(res$alpha, 0.4)
})

test_that("tasne", {
  res <- smallvis(iris10, Y_init = iris10_Y, method = "tasne", eta = 0.1,
                  perplexity = 5,
                  epoch_callback = NULL, verbose = FALSE, ret_extra = TRUE)
  expect_equal(res$Y, c2y(-9.398, 6.729, 5.751, 9.675, -8.758, -23.45, 4.814, -4.545,
                        16.18, 2.997, 2.331, 4.598, -1.154, -0.7406, -2.752, -0.8487,
                        -5.518, 1.061, -0.7827, 3.806), tolerance = 1e-4)
  expect_equal(final_cost(res), 0.5719, tolerance = 1e-4)
})

test_that("jse", {
  res <- smallvis(iris10, Y_init = iris10_Y, method = "jse", eta = 0.1,
                  perplexity = 5,
                  epoch_callback = NULL, verbose = FALSE, ret_extra = TRUE)
  expect_equal(res$Y, c2y(-0.8728, 0.3486, 0.6272, 0.9292, -0.9257, -2.032, 0.5479, -0.5908,
                        1.783, 0.185, -0.4006, 0.6197, -0.2429, -0.06082, -0.5231, 0.8692,
                        -0.848, -0.2266, 0.2949, 0.5181), tolerance = 1e-4)
  expect_equal(final_cost(res), 0.438, tolerance = 1e-5)
})

test_that("nerv", {
  res <- smallvis(iris10, Y_init = iris10_Y, method = "nerv", eta = 0.1,
                  perplexity = 5,
                  epoch_callback = NULL, verbose = FALSE, ret_extra = TRUE)
  expect_equal(res$Y, c2y(-0.8838, 0.3817, 0.5606, 0.8724, -0.9276, -1.791, 0.4291, -0.5919,
                        1.732, 0.2189, -0.3627, 0.6106, -0.2454, -0.09668, -0.4693, 0.927,
                        -0.7941, -0.2128, 0.1265, 0.5168), tolerance = 1e-4)
  expect_equal(final_cost(res), 0.5399, tolerance = 1e-5)
})

test_that("ee", {
  res <- smallvis(iris10, Y_init = iris10_Y, method = "ee", eta = 0.1,
                  perplexity = 5,
                  epoch_callback = NULL, verbose = FALSE, ret_extra = TRUE)
  expect_equal(res$Y, c2y(-1.785, 1.447, 0.8011, 1.907, -1.644, -3.393, 0.2531, -0.8798,
                        3.017, 0.276, 0.8166, 1.267, -0.2262, -0.8756, -0.9315, -0.1936,
                        -1.593, 0.1561, 0.1982, 1.383), tolerance = 1e-3)
  expect_equal(final_cost(res), 5.563, tolerance = 1e-4)
})

test_that("umap", {
  res <- smallvis(iris10, Y_init = iris10_Y, method = "umap", eta = 0.1,
                  perplexity = 5,
                  epoch_callback = NULL, verbose = FALSE, ret_extra = TRUE)
  expect_equal(res$Y, c2y(-1.297, 0.5526, 0.8944, 1.112, -1.337, -1.578, 0.824, -1.001,
                        1.41, 0.4207, -0.006171, 0.3732, -0.2137, -0.05575, -0.02374,
                        -0.03407, -0.4788, 0.02117, 0.05909, 0.3587), tolerance = 1e-3)
  expect_equal(final_cost(res), 6.144, tolerance = 1e-4)
})

test_that("largevis", {
  res <- smallvis(iris10, Y_init = iris10_Y, method = "largevis", eta = 0.1,
                  perplexity = 5,
                  epoch_callback = NULL, verbose = FALSE, ret_extra = TRUE)
  expect_equal(res$Y, c2y(-8.834, 5.819, 3.866, 7.954, -7.697, -13.68, 2.955, -4.704,
                        12.38, 1.944, 1.879, 4.428, -1.165, -1.434, -2.13, -1.525, -5.295,
                        1.379, -0.3885, 4.251), tolerance = 1e-4)
  expect_equal(final_cost(res), 46.99, tolerance = 1e-4)
})

