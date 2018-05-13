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

  expected_itercosts <- c(0.3087, 0.04999, 0.04132, 0.04129, 0.04129)
  names(expected_itercosts) <- c(100, 200, 300, 400, 500)
  expect_equal(i10_tsne$itercosts, expected_itercosts, tolerance = 1e-4)

  expect_equal(i10_tsne$costs, c(0.005045, 0.0005255, -0.005447, 0.003863,
                                 0.0005298,  0.006410,  0.01191, 0.007799,
                                 0.003098,  0.007556), tolerance = 1e-5)

  expect_equal(i10_tsne$Y, expected_Y, tolerance = 0.1)
})

test_that("extra return values and iter0 cost", {
  i10_tsne <- smallvis(iris10, Y_init = iris10_Y, perplexity = 5,
                       epoch_callback = NULL, verbose = FALSE,
                       ret_extra = TRUE, iter0_cost = TRUE)
  expect_equal(i10_tsne$N, 10)
  expect_equal(i10_tsne$origD, 4)
  expect_equal(i10_tsne$scale, "absmax")
  expect_equal(i10_tsne$Y_init, "matrix")
  expect_equal(i10_tsne$method, "tsne")
  expect_equal(i10_tsne$perplexity, 5)

  expected_itercosts <- c(0.3542, 0.3087, 0.04999, 0.04132, 0.04129, 0.04129)
  names(expected_itercosts) <- c(0, 100, 200, 300, 400, 500)
  expect_equal(i10_tsne$itercosts, expected_itercosts, tolerance = 1e-4)

  expect_equal(i10_tsne$costs, c(0.005045, 0.0005255, -0.005447, 0.003863,
                                 0.0005298,  0.006410,  0.01191, 0.007799,
                                 0.003098,  0.007556), tolerance = 1e-5)

  expect_equal(i10_tsne$Y, expected_Y, tolerance = 0.1)
})

test_that("early and late exaggeration", {
  # P should be a multiple when exaggeration is on
  i10_tsne <- smallvis(iris10, Y_init = iris10_Y, perplexity = 5,
                       epoch_callback = NULL, verbose = FALSE,
                       ret_extra = c("P"), exaggeration_factor = 4,
                       stop_lying_iter = 11, max_iter = 10)
  expect_equal(sum(i10_tsne$P), 4)
  expect_equal(i10_tsne$exaggeration_factor, 4)
  expect_equal(i10_tsne$stop_lying_iter, 11)
  # Don't report late exaggeration if not used
  expect_null(i10_tsne$start_late_lying_iter)
  expect_null(i10_tsne$late_exaggeration_factor)

  # P should be back to normal when exaggeration is off
  i10_tsne <- smallvis(iris10, Y_init = iris10_Y, perplexity = 5,
                       epoch_callback = NULL, verbose = FALSE,
                       ret_extra = c("P"), exaggeration_factor = 4,
                       stop_lying_iter = 5, max_iter = 10)
  expect_equal(sum(i10_tsne$P), 1)
  expect_equal(i10_tsne$exaggeration_factor, 4)
  expect_equal(i10_tsne$stop_lying_iter, 5)
  expect_null(i10_tsne$start_late_lying_iter)
  expect_null(i10_tsne$late_exaggeration_factor)


  # P should be back to a multiple with late exaggeration
  i10_tsne <- smallvis(iris10, Y_init = iris10_Y, perplexity = 5,
                       epoch_callback = NULL, verbose = FALSE,
                       ret_extra = c("P"),
                       exaggeration_factor = 4, stop_lying_iter = 5,
                       late_exaggeration_factor = 4, start_late_lying_iter = 9,
                       max_iter = 10)
  expect_equal(sum(i10_tsne$P), 4)
  expect_equal(i10_tsne$exaggeration_factor, 4)
  expect_equal(i10_tsne$late_exaggeration_factor, 4)
  expect_equal(i10_tsne$stop_lying_iter, 5)
  expect_equal(i10_tsne$start_late_lying_iter, 9)

  # different early and late exaggeration
  i10_tsne <- smallvis(iris10, Y_init = iris10_Y, perplexity = 5,
                       epoch_callback = NULL, verbose = FALSE,
                       ret_extra = c("P"), exaggeration_factor = 4,
                       late_exaggeration_factor = 2,
                       stop_lying_iter = 5,  start_late_lying_iter = 9,
                       max_iter = 10)
  expect_equal(sum(i10_tsne$P), 2)
  expect_equal(i10_tsne$exaggeration_factor, 4)
  expect_equal(i10_tsne$late_exaggeration_factor, 2)
  expect_equal(i10_tsne$stop_lying_iter, 5)
  expect_equal(i10_tsne$start_late_lying_iter, 9)
})


test_that("mmds", {
  res <- smallvis(iris10, Y_init = iris10_Y, method = "mmds", eta = 0.1,
                       epoch_callback = NULL, verbose = FALSE, ret_extra = TRUE)
  expect_equal(res$Y, c2y(0.5031, -0.3515, -0.3407, -0.59, 0.4913, 1.463, -0.209, 0.2736,
                            -1.056, -0.1838, -0.02805, -0.4261, 0.1157, 0.01604, 0.2117,
                            0.0396, 0.4632, -0.09042, 0.0585, -0.3602), tolerance = 1e-4)
  expect_equal(final_cost(res), 0.1715, tolerance = 1e-5)
})

test_that("mmds without scaling", {
  res <- smallvis(iris10, Y_init = iris10_Y, method = "mmds", scale = FALSE,
                  eta = 0.1,
                  epoch_callback = NULL, verbose = FALSE, ret_extra = TRUE)
  expect_equal(res$Y, c2y(-0.2962, 0.2143, 0.1991, 0.3477, -0.2932, -0.8634, 0.1157,
                          -0.1599, 0.6217, 0.1143, 0.02476, 0.2456, -0.07381, -0.01909,
                          -0.1169, 0.0005308, -0.2766, 0.0578, -0.05174, 0.2095), tolerance = 1e-4)
  expect_equal(final_cost(res), 0.0597, tolerance = 1e-5)
})

test_that("mmds with distance matrix", {
  iris10d <- dist(iris10)
  res <- smallvis(iris10d, Y_init = iris10_Y, method = "mmds", scale = FALSE,
                  eta = 0.1,
                  epoch_callback = NULL, verbose = FALSE, ret_extra = TRUE)
  expect_equal(res$Y, c2y(-0.2962, 0.2143, 0.1991, 0.3477, -0.2932, -0.8634, 0.1157,
                          -0.1599, 0.6217, 0.1143, 0.02476, 0.2456, -0.07381, -0.01909,
                          -0.1169, 0.0005308, -0.2766, 0.0578, -0.05174, 0.2095), tolerance = 1e-4)
  expect_equal(final_cost(res), 0.0597, tolerance = 1e-5)
})

test_that("sammon", {
  res <- smallvis(iris10, Y_init = iris10_Y, method = "sammon", eta = 0.1,
                  epoch_callback = NULL, verbose = FALSE, ret_extra = TRUE)
  expect_equal(res$Y, c2y(-0.5143, 0.3886, 0.3202, 0.6053, -0.4982, -1.473, 0.2011, -0.2667,
                            1.061, 0.1761, 0.02941, 0.4171, -0.1152, -0.0331, -0.2058, 0.03922,
                            -0.4762, 0.09549, -0.1027, 0.3518), tolerance = 1e-3)
  expect_equal(final_cost(res), 0.003427, tolerance = 1e-5)
})

test_that("smmds", {
  res <- smallvis(iris10, Y_init = iris10_Y, method = "smmds", eta = 0.001,
                  epoch_callback = NULL, verbose = FALSE, ret_extra = TRUE)
  expect_equal(res$Y, c2y(-0.4952, 0.3507, 0.3716, 0.5528, -0.502, -1.449, 0.1852, -0.2768,
                            1.042, 0.2206, 0.1035, 0.4062, -0.1173, -0.07106, -0.1869, -0.02864,
                            -0.459, 0.09674, -0.1043, 0.3608), tolerance = 1e-3)
  expect_equal(final_cost(res), 0.3038, tolerance = 1e-4)
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

test_that("ballmmds", {
  res <- smallvis(iris10, Y_init = iris10_Y, method = list("ballmmds", f = 0.5), eta = 0.1,
                  epoch_callback = NULL, verbose = FALSE, ret_extra = TRUE, max_iter = 50)
  expect_equal(res$Y, c2y(0.0849, 1.09, 0.9095, 1.208, 0.05264, -6.899, 0.7066, 0.3632,
                          1.624, 0.8593, -0.1728, -0.381, 0.1399, 0.1256, 0.05778, -0.1192,
                          0.4718, -0.1738, 0.3942, -0.3424), tolerance = 1e-4)
  expect_equal(final_cost(res), 0.08589, tolerance = 1e-4)
})


test_that("knnmmds", {
  res <- smallvis(iris10, Y_init = iris10_Y, method = "knnmmds", eta = 0.1, perplexity = 5,
                  epoch_callback = NULL, verbose = FALSE, ret_extra = TRUE, max_iter = 50)
  expect_equal(res$Y, c2y(-0.1916, -0.2973, 0.1988, 0.228, 0.03356, -0.419, 0.4968, -0.201,
                        0.4247, -0.2729, -0.4865, 0.5189, 0.285, 0.5853, -0.5463, -1.44,
                        0.03945, -0.2379, 0.9893, 0.2931), tolerance = 1e-3)
  expect_equal(final_cost(res), 0.1155, tolerance = 1e-4)
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
  expect_equal(res$Y, c2y(-1.359, 0.7321, 0.8316, 1.071, -1.346, -1.606, 0.7614, -1.036,
                          1.407, 0.543, -0.0262, 0.5369, -0.1322, -0.1236, -0.07116, -0.08141,
                          -0.4839, 0.01953, -0.1767, 0.5387), tolerance = 1e-3)
  expect_equal(final_cost(res), 13.15, tolerance = 1e-4)
})

test_that("tumap", {
  res <- smallvis(iris10, Y_init = iris10_Y, method = "tumap", eta = 0.1,
                  perplexity = 5, max_iter = 100,
                  epoch_callback = NULL, verbose = FALSE, ret_extra = TRUE)
  expect_equal(res$Y, c2y(-1.609, 0.8876, 0.9294, 1.305, -1.586, -1.982, 0.7682, -1.107,
                          1.792, 0.6012, -0.002293, 0.7693, -0.2082, -0.2003, -0.0997,
                          -0.1114, -0.701, 0.0876, -0.2937, 0.7598), tolerance = 1e-3)
  expect_equal(final_cost(res), 11.28, tolerance = 1e-4)
})

test_that("ntumap", {
  res <- smallvis(iris10, Y_init = iris10_Y, method = "ntumap", eta = 10,
                  perplexity = 5,
                  epoch_callback = NULL, verbose = FALSE, ret_extra = TRUE)
  expect_equal(res$Y, c2y(-14.05, 10.35, 5.266, 11.03, -11.35, -18.74, 3.344, -6.781,
                          17.74, 3.193, 3.335, 6.7, -0.7754, -3.755, -3.275, -2.134, -7.851,
                          2.971, -2.838, 7.622), tolerance = 1e-3)
  expect_equal(final_cost(res), 0.1622, tolerance = 1e-4)
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

test_that("tsne with L-BFGS", {
  res <- smallvis::smallvis(iris10, Y_init = iris10_Y, perplexity = 5,
                       epoch_callback = NULL, verbose = FALSE,
                       opt = list("L-BFGS"), ret_extra = TRUE, max_iter = 25)
  expect_equal(res$Y, c2y(-3.534, 2.266, 1.774, 3.28, -3.286, -5.684, 1.287, -2.035,
                           5.02, 0.9119, 0.7444, 1.769, -0.6249, -0.4774, -0.7323, -0.2723,
                           -2.118, 0.4957, -0.3813, 1.598), tolerance = 1e-3)
  expect_equal(final_cost(res), 0.02590, tolerance = 1e-4)
})

test_that("Miscellany", {
  res <- smallvis(iris10, Y_init = iris10_Y, method = "bnerv", eta = 0.1,
                  perplexity = 5,
                  epoch_callback = NULL, verbose = FALSE, ret_extra = TRUE)
  expect_equal(res$Y, c2y(-0.4777, 0.3406, 0.3258, 0.574, -0.4849, -1.468, 0.2147, -0.2627,
                          1.012, 0.2258, 0.02774, 0.3532, -0.1, -0.05147, -0.1651, 0.06006,
                          -0.4242, 0.07663, -0.09239, 0.3156), tolerance = 1e-3)
  expect_equal(final_cost(res), 0.0574, tolerance = 1e-4)


  testthat::expect_true(!is.null(res))
  res <- smallvis(iris10, Y_init = iris10_Y, method = "rsrnerv", eta = 0.1,
                  perplexity = 5,
                  epoch_callback = NULL, verbose = FALSE, ret_extra = TRUE)
  expect_equal(res$Y, c2y(-0.9464, 0.4992, 0.4742, 0.8761, -0.9257, -1.231, 0.2172, -0.6745,
                          1.411, 0.2996, 0.1939, 0.662, -0.2914, -0.177, 0.02104, -0.5248,
                          -0.7185, 0.2888, -0.09876, 0.6447), tolerance = 1e-3)
  expect_equal(final_cost(res), 0.1731, tolerance = 1e-4)

  res <- smallvis(iris10, Y_init = iris10_Y, method = "rsrjse", eta = 0.1,
                  perplexity = 5,
                  epoch_callback = NULL, verbose = FALSE, ret_extra = TRUE)
  expect_equal(res$Y, c2y(-0.9463, 0.5184, 0.4731, 0.884, -0.9286, -1.281, 0.2148, -0.6778,
                          1.429, 0.3148, 0.1945, 0.6732, -0.2905, -0.1792, 0.02602, -0.5414,
                          -0.7223, 0.2907, -0.1077, 0.6568), tolerance = 1e-3)
  expect_equal(final_cost(res), 0.1627, tolerance = 1e-4)

  res <- smallvis(iris10, Y_init = iris10_Y, method = "btsne", eta = 0.1,
                  perplexity = 5,
                  epoch_callback = NULL, verbose = FALSE, ret_extra = TRUE)
  expect_equal(res$Y, c2y(-1.437, 0.949, 0.8898, 1.374, -1.401, -2.975, 0.7103, -0.813,
                          2.22, 0.4825, 0.3022, 0.7113, -0.2084, -0.1298, -0.3676, -0.06477,
                          -0.7179, 0.1204, -0.2042, 0.5588), tolerance = 1e-3)
  expect_equal(final_cost(res), 0.07386, tolerance = 1e-4)

  res <- smallvis(iris10, Y_init = iris10_Y, method = "bssne", eta = 0.1,
                  perplexity = 5,
                  epoch_callback = NULL, verbose = FALSE, ret_extra = TRUE)
  expect_equal(res$Y, c2y(-0.4612, 0.2633, 0.2751, 0.4448, -0.4517, -0.9009, 0.1995,
                          -0.2956, 0.7688, 0.1578, -0.07989, 0.3133, -0.08017, -0.03674,
                          -0.1972, 0.1591, -0.3072, 0.01725, -0.06213, 0.2736), tolerance = 1e-3)
  expect_equal(final_cost(res), 0.1089, tolerance = 1e-4)

  res <- smallvis(iris10, Y_init = iris10_Y, method = "basne", eta = 0.1,
                  perplexity = 5,
                  epoch_callback = NULL, verbose = FALSE, ret_extra = TRUE)
  expect_equal(res$Y, c2y(-0.4776, 0.3416, 0.3254, 0.5736, -0.4853, -1.467, 0.2133, -0.2623,
                          1.012, 0.2268, 0.02933, 0.352, -0.1008, -0.05319, -0.1635, 0.06401,
                          -0.4247, 0.07757, -0.09564, 0.3148), tolerance = 1e-3)
  expect_equal(final_cost(res), 2.123, tolerance = 1e-4)

  res <- smallvis(iris10, Y_init = iris10_Y, method = "btasne", eta = 0.1,
                  perplexity = 5,
                  epoch_callback = NULL, verbose = FALSE, ret_extra = TRUE)
  expect_equal(res$Y, c2y(-6.105, 4.355, 3.715, 6.291, -5.661, -15.21, 3.147, -2.944,
                          10.5, 1.911, 1.465, 3.026, -0.7081, -0.4368, -1.832, -0.6734,
                          -3.566, 0.6635, -0.4344, 2.497), tolerance = 1e-3)
  expect_equal(final_cost(res), 104.7, tolerance = 1e-4)

  res <- smallvis(iris10, Y_init = iris10_Y, method = "trmsne", eta = 0.1,
                  perplexity = 5,
                  epoch_callback = NULL, verbose = FALSE, ret_extra = TRUE)
  expect_equal(res$Y, c2y(-3.487, 2.223, 1.753, 3.228, -3.229, -5.592, 1.276, -2.013,
                           4.941, 0.8986, 0.7316, 1.745, -0.6139, -0.4619, -0.7119, -0.3068,
                           -2.083, 0.504, -0.3769, 1.574), tolerance = 1e-3)
  expect_equal(final_cost(res), 0.1207, tolerance = 1e-4)

  res <- smallvis(iris10, Y_init = iris10_Y, method = "tmsne", eta = 0.1,
                  perplexity = 5,
                  epoch_callback = NULL, verbose = FALSE, ret_extra = TRUE)
  expect_equal(res$Y, c2y(-2.699, 1.629, 1.619, 3.104, -2.543, -5.294, 0.8766, -1.73,
                          4.422, 0.6155, 0.4124, 1.596, -0.6119, 0.1443, -0.4072, -0.2042,
                          -2.248, 0.5002, -0.6253, 1.444), tolerance = 1e-3)
  expect_equal(final_cost(res), 0.1226, tolerance = 1e-4)

  res <- smallvis(iris10, Y_init = iris10_Y, method = "trsrsne", eta = 0.1,
                  perplexity = 5,
                  epoch_callback = NULL, verbose = FALSE, ret_extra = TRUE)
  expect_equal(res$Y, c2y(-3.524, 2.286, 1.786, 3.237, -3.254, -5.546, 1.27, -2.062,
                          4.891, 0.9163, 0.6728, 1.785, -0.5903, -0.4569, -0.7359, -0.3572,
                          -2.115, 0.4887, -0.3448, 1.654), tolerance = 1e-3)
  expect_equal(final_cost(res), 0.2683, tolerance = 1e-4)

  res <- smallvis(iris10, Y_init = iris10_Y, method = "arsrsne", eta = 0.1,
                  perplexity = 5,
                  epoch_callback = NULL, verbose = FALSE, ret_extra = TRUE)
  expect_equal(res$Y, c2y(-0.9463, 0.4984, 0.4736, 0.8745, -0.925, -1.225, 0.2149, -0.6731,
                          1.408, 0.2993, 0.1963, 0.659, -0.2926, -0.1789, 0.02274, -0.5197,
                          -0.7179, 0.2904, -0.1017, 0.6423), tolerance = 1e-3)
  expect_equal(final_cost(res), 0.1741, tolerance = 1e-4)
})

test_that("repeated runs", {
  # Random intialization so don't care about numerical results
  # Test that combinations of keep_all and ret_extra return correct types

  res_keep_ret <- smallvis_rep(
    n = 3, keep_all = TRUE, X = iris10, scale = FALSE, perplexity = 5,
    ret_extra = TRUE, max_iter = 5, epoch_callback = NULL, verbose = FALSE)
  expect_length(res_keep_ret, 3)
  expect_is(res_keep_ret[[1]], "list")
  expect_length(res_keep_ret[[1]]$all_costs, 3)
  expect_equal(res_keep_ret[[1]]$best_rep, which.min(res_keep_ret[[1]]$all_costs))

  res_nokeep_ret <- smallvis_rep(
    n = 3, keep_all = FALSE, X = iris10, scale = FALSE, perplexity = 5,
    ret_extra = TRUE, max_iter = 5, epoch_callback = NULL, verbose = FALSE)
  expect_is(res_nokeep_ret, "list")
  expect_length(res_nokeep_ret$all_costs, 3)
  expect_equivalent(res_nokeep_ret$itercosts[length(res_nokeep_ret$itercosts)],
               res_nokeep_ret$all_costs[which.min(res_nokeep_ret$all_costs)])

  res_keep_noret <- smallvis_rep(
    n = 3, keep_all = TRUE, X = iris10, scale = FALSE, perplexity = 5,
    ret_extra = FALSE, max_iter = 5, epoch_callback = NULL, verbose = FALSE)
  expect_is(res_keep_noret, "list")
  expect_length(res_keep_noret, 3)
  expect_is(res_keep_noret[[1]], "matrix")

  res_nokeep_noret <- smallvis_rep(
    n = 3, keep_all = FALSE, X = iris10, scale = FALSE, perplexity = 5,
    ret_extra = FALSE, max_iter = 5, epoch_callback = NULL, verbose = FALSE)
  expect_is(res_nokeep_noret, "matrix")
})

test_that("multiple equal perplexities same as one perplexity", {
  i10_tsne <- smallvis(iris10, Y_init = iris10_Y, perplexity = rep(5, 10),
                       epoch_callback = NULL, verbose = FALSE)
  expect_equal(i10_tsne, expected_Y, tolerance = 0.1)
})

test_that("multiple perplexities", {
  i10_tsne <- smallvis(iris10, Y_init = iris10_Y, perplexity = c(rep(5, 5), rep(6, 5)),
                       epoch_callback = NULL, verbose = FALSE)
  expect_equal(i10_tsne, c2y(-335.8, 169.9, 354.1, 250, -412, -629.5, 511.4, -193.4, 237.1,
                           48.18, 271.5, 120.7, -35.23, -208.6, 87.18, 284, -171.8, 121.4,
                           -436.3, -32.79), tolerance = 1e-3)
})
