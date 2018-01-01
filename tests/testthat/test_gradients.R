library(smallvis)
context("Gradients")

perp <- 5

test_that("SNE", {
  test_grad("tsne", perplexity = perp)
  test_grad("ssne", perplexity = perp)
  test_grad("asne", perplexity = perp)
  test_grad("wtsne", perplexity = perp)
})

test_that("HSSNE", {
  # fd starts losing accuracy for alpha = 0
  test_grad("hssne", perplexity = perp, alpha = 1e-6)
  test_grad("hssne", perplexity = perp, alpha = 0.5)
  test_grad("hssne", perplexity = perp, alpha = 1)
})

test_that("JSE", {
  test_grad("jse", kappa = 0.5, perplexity = perp)
  # fd starts losing accuracy for kappa = 0
  test_grad("jse", kappa = 1e-5, perplexity = perp)
  test_grad("jse", kappa = 1, perplexity = perp)
})

test_that("NeRV", {
  test_grad("nerv", lambda = 0.5, perplexity = perp)
  test_grad("nerv", lambda = 0, perplexity = perp)
  test_grad("nerv", lambda = 1, perplexity = perp)
})

test_that("MDS", {
  test_grad("mmds")
  test_grad("smmds")
  test_grad("sammon")
  test_grad("geommds", k = 3)
})

test_that("EE", {
  test_grad("ee", lambda = 1, perplexity = perp)
  test_grad("ee", lambda = 100, perplexity = perp)
  test_grad("ee", lambda = 1000, perplexity = perp)
})

test_that("LargeVis", {
  # LV gradients can be extremely large compared to other methods
  test_grad("largevis", tolerance = 1e-3, gamma = 7, perplexity = perp, gr_eps = 0)
  test_grad("largevis", gamma = 1e-3, perplexity = perp, gr_eps = 0)
  test_grad("largevis", tolerance = 1e-2, gamma = 1000, perplexity = perp, gr_eps = 0)
})

test_that("UMAP", {
  test_grad("umap", perplexity = perp, gr_eps = 0)
  test_grad("umap", tolerance = 1e-5, spread = 10, min_dist = 0.01, perplexity = perp, gr_eps = 0)
  test_grad("umap", spread = 0.5, min_dist = 1e-4, perplexity = perp, gr_eps = 0)
  test_grad("tumap", perplexity = perp, gr_eps = 0)
  test_grad("ntumap", perplexity = perp, gr_eps = 0)
})

