library(smallvis)
context("Intrinsic Dimensionality Perplexity")

test_that("Default trial perplexities are size appropriate", {
expect_equal(idp_perps(nrow(iris10)), c(4))
expect_equal(idp_perps(nrow(iris)), c(4, 8, 16, 32, 64))
expect_equal(idp_perps(6000), c(4, 8, 16, 32, 64, 128))
})

test_that("Can use IDP with one default perplexity", {
  i10_tsne_idp <- smallvis(iris10, Y_init = iris10_Y, perplexity = "idp",
                       epoch_callback = NULL, verbose = FALSE,
                       ret_extra = TRUE)
  expect_equal(i10_tsne_idp$perplexity, 4)
})

test_that("Can provide perplexities to IDP", {
  i10_tsne_idp <- smallvis(iris10, Y_init = iris10_Y,
                           perplexity = list("idp", 2:4),
                           epoch_callback = NULL, verbose = FALSE,
                           ret_extra = TRUE)
  expect_equal(i10_tsne_idp$perplexity, 3)
})
