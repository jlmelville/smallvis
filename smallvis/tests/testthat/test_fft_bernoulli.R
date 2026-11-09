test_that("FFT Bernoulli methods resolve short forces on a wide grid", {
  centers <- as.matrix(expand.grid(x = 5 * (1:8), y = 5 * (1:8)))
  n <- nrow(centers)
  angle <- seq_len(n)
  distance <- exp(seq(log(1e-4), log(0.8), length.out = n))
  Y <- unname(rbind(
    centers,
    centers + cbind(distance * cos(angle), distance * sin(angle)),
    c(0, 0),
    c(100, 100)
  ))
  prepare <- function(method, ...) {
    cost <- create_cost(
      method,
      perplexity = 4,
      eps = 1e-9,
      n_threads = 2,
      use_cpp = FALSE,
      theta = 0,
      inp_kernel = NULL,
      nn = NULL,
      bh = FALSE
    )
    # Isolate repulsion from the sparse attraction in this fixture.
    cost$P <- Matrix::sparseMatrix(
      i = integer(),
      j = integer(),
      x = numeric(),
      dims = rep(nrow(Y), 2L)
    )
    if (method[[1L]] == "umap") {
      ab <- find_ab_params()
      cost$a <- ab[["a"]]
      cost$b <- ab[["b"]]
    }
    cost_cache_input(cost)
  }
  for (name in c("largevis", "tumap", "umap")) {
    for (stabilization in c(0.1, 0)) {
      method <- list(name, gr_eps = stabilization)
      exact <- prepare(c(method, list(repulsion = "bh")))
      fft <- prepare(c(method, list(repulsion = "fft")))
      G <- cost_grad(exact, Y)$G
      approx <- cost_grad(fft, Y)$G
      expect_lt(sqrt(sum((G - approx)^2) / sum(G^2)), 0.03)
      expect_equal(
        cost_point(fft, Y)$pcost,
        cost_point(exact, Y)$pcost,
        tolerance = 0.01
      )
    }
  }
})

test_that("FFT Bernoulli methods share the saved graph and preserve collisions", {
  X <- cbind(seq_len(24), sin(seq_len(24)), cos(seq_len(24)))
  neighbors <- get_nn(stats::dist(X), k = 13L)
  Y <- cbind(sin(seq_len(24)), cos(seq_len(24)))
  Y[2L, ] <- Y[1L, ]
  local_mocked_bindings(get_nn = function(...) stop("neighbor search ran"))
  prepare <- function(method, ...) {
    smallvis_run(
      X,
      method = method,
      neighbors = neighbors,
      perplexity = 4,
      scale = FALSE,
      prepare_only = TRUE,
      epoch_callback = FALSE,
      verbose = FALSE,
      ...
    )$cost_fn
  }
  for (name in c("largevis", "tumap", "umap")) {
    fft <- prepare(list(name, repulsion = "fft"))
    bh <- prepare(name, bh = TRUE, theta = 0)
    expect_identical(fft$P, bh$P)
    expect_equal(cost_grad(fft, Y)$G, cost_grad(bh, Y)$G, tolerance = 0.01)
    compact <- Y / 100
    expect_equal(
      cost_grad(fft, compact)$G,
      cost_grad(bh, compact)$G,
      tolerance = 0.01
    )
    expect_equal(
      cost_point(fft, Y)$pcost,
      cost_point(bh, Y)$pcost,
      tolerance = 0.01
    )
  }
  unregularized <- prepare(list("largevis", repulsion = "fft", gr_eps = 0))
  expect_error(
    cost_grad(unregularized, Y)$G,
    "undefined repulsion at zero distance"
  )
})
