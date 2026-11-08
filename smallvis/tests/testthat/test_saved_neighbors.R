saved_neighbors_fixture <- function() {
  X <- cbind(seq_len(24), sin(seq_len(24)), cos(seq_len(24)))
  D <- stats::dist(X)
  list(X = X, D = D, graph = get_nn(D, k = 10L))
}

test_that("dense and Barnes-Hut methods can use the supplied graph without search", {
  fixture <- saved_neighbors_fixture()
  Y <- cbind(sin(seq_len(24)), cos(seq_len(24))) / 3
  run <- function(method, kernel, bh, X, ...) {
    smallvis(
      X,
      method = method,
      bh = bh,
      inp_kernel = kernel,
      perplexity = 3,
      Y_init = Y,
      theta = 0,
      max_iter = 1,
      ret_extra = c("P", "beta"),
      epoch_callback = FALSE,
      verbose = FALSE,
      ...
    )
  }
  cases <- expand.grid(
    method = c("tsne", "hssne", "tee", "largevis", "tumap", "jssne"),
    kernel = c("knn", "nngaussian"),
    stringsAsFactors = FALSE
  )
  cases <- rbind(
    transform(cases, bh = TRUE),
    data.frame(method = c("tsne", "umap", "tumap"), kernel = "skd", bh = TRUE),
    data.frame(
      method = c("ssne", "hlsne", "btsne", "umap"),
      kernel = c("knn", "nngaussian", "nngaussian", "skd"),
      bh = FALSE
    )
  )
  expected <- lapply(seq_len(nrow(cases)), function(i) {
    run(cases$method[i], cases$kernel[i], cases$bh[i], fixture$D, nn = "exact")
  })
  local_mocked_bindings(get_nn = function(...) stop("neighbor search ran"))
  for (i in seq_len(nrow(cases))) {
    graph <- fixture$graph
    # Test both omitted self and self outside the first column.
    cols <- switch(cases$kernel[i], knn = 2:10, skd = 2:3, c(2:10, 1L))
    graph$idx <- graph$idx[, cols]
    graph$dist <- graph$dist[, cols]
    if (cases$method[i] == "tee") {
      storage.mode(graph$idx) <- "double"
    }
    actual <- run(
      cases$method[i],
      cases$kernel[i],
      cases$bh[i],
      fixture$X,
      neighbors = graph
    )
    expect_equal(actual$P, expected[[i]]$P)
    expect_equal(actual$Y, expected[[i]]$Y)
    expect_equal(actual$costs, expected[[i]]$costs)
    if (cases$method[i] == "btsne") {
      expect_length(actual$beta, nrow(fixture$X))
      expect_equal(actual$beta, expected[[i]]$beta)
    }
  }
})

test_that("exact Barnes-Hut LargeVis matches dense forces on the saved graph", {
  fixture <- saved_neighbors_fixture()
  Y <- cbind(sin(seq_len(24)), cos(seq_len(24))) / 3
  for (gr_eps in c(0, 0.2)) {
    method <- list(
      "largevis",
      gamma = 0.03,
      gr_eps = gr_eps,
      normalize = FALSE,
      row_weight = seq_len(24) / 24
    )
    coordinates <- Y
    if (gr_eps != 0) {
      coordinates[2, ] <- coordinates[1, ]
    }
    local({
      prepare <- function(X, bh, ...) {
        smallvis_run(
          X,
          method = method,
          bh = bh,
          inp_kernel = if (gr_eps == 0) "knn" else "nngaussian",
          perplexity = 3,
          theta = 0,
          eps = if (gr_eps == 0) 0 else 1e-5,
          n_threads = 2,
          prepare_only = TRUE,
          epoch_callback = FALSE,
          verbose = FALSE,
          ...
        )$cost_fn
      }
      dense <- prepare(fixture$D, FALSE)
      local_mocked_bindings(get_nn = function(...) stop("neighbor search ran"))
      bh <- prepare(fixture$X, TRUE, neighbors = fixture$graph)
      expect_equal(as.matrix(bh$P), dense$P, tolerance = 1e-12)
      expect_equal(
        cost_grad(bh, coordinates)$G,
        cost_grad(dense, coordinates)$G,
        tolerance = 1e-12
      )
      expect_equal(
        cost_point(bh, coordinates)$pcost,
        sum(cost_point(dense, coordinates)$pcost),
        tolerance = 1e-12
      )
      if (gr_eps == 0) {
        coordinates[2, ] <- coordinates[1, ]
        expect_error(
          cost_grad(bh, coordinates),
          "undefined repulsion at zero distance"
        )
      }
    })
  }
})

test_that("exact Barnes-Hut t-UMAP matches the dense objective on the saved graph", {
  fixture <- saved_neighbors_fixture()
  Y <- cbind(sin(seq_len(24)), cos(seq_len(24))) / 3
  for (kernel in c("knn", "nngaussian", "skd")) {
    method <- if (kernel == "knn") {
      list("tumap", gr_eps = 0)
    } else {
      list("tumap", gr_eps = 0.2, row_weight = seq_len(24) / 24)
    }
    coordinates <- Y
    if (method$gr_eps != 0) {
      coordinates[2, ] <- coordinates[1, ]
    }
    local({
      prepare <- function(X, bh, ...) {
        smallvis_run(
          X,
          method = method,
          bh = bh,
          inp_kernel = kernel,
          perplexity = if (kernel == "skd") 4.5 else 3,
          theta = 0,
          eps = 1e-5,
          n_threads = 2,
          prepare_only = TRUE,
          epoch_callback = FALSE,
          verbose = FALSE,
          ...
        )$cost_fn
      }
      dense <- prepare(fixture$D, FALSE)
      local_mocked_bindings(get_nn = function(...) stop("neighbor search ran"))
      bh <- prepare(fixture$X, TRUE, neighbors = fixture$graph)
      expect_s4_class(bh$P, "dgCMatrix")
      expect_equal(as.matrix(bh$P), dense$P, tolerance = 1e-12)
      expect_identical(bh$gradient_only, dense$gradient_only)
      expect_equal(
        cost_grad(bh, coordinates)$G,
        cost_grad(dense, coordinates)$G,
        tolerance = 1e-12
      )
      expect_equal(
        cost_point(bh, coordinates)$pcost,
        sum(cost_point(dense, coordinates)$pcost),
        tolerance = 1e-12
      )
    })
  }
})

test_that("exact Barnes-Hut UMAP matches dense cost and gradient on a saved graph", {
  fixture <- saved_neighbors_fixture()
  Y <- cbind(sin(seq_len(24)), cos(seq_len(24))) / 3
  for (min_dist in c(0.001, 0.5)) {
    method <- list("umap", min_dist = min_dist, gr_eps = 0)
    local({
      prepare <- function(X, bh, ...) {
        smallvis_run(
          X,
          method = method,
          bh = bh,
          perplexity = 4.5,
          theta = 0,
          eps = 1e-5,
          prepare_only = TRUE,
          epoch_callback = FALSE,
          verbose = FALSE,
          ...
        )$cost_fn
      }
      dense <- prepare(fixture$D, FALSE)
      local_mocked_bindings(get_nn = function(...) stop("neighbor search ran"))
      bh <- prepare(fixture$X, TRUE, neighbors = fixture$graph)
      expect_identical(bh$a, dense$a)
      expect_identical(bh$b, dense$b)
      expect_equal(as.matrix(bh$P), dense$P, tolerance = 1e-12)
      expect_equal(
        cost_grad(bh, Y)$G,
        cost_grad(dense, Y)$G,
        tolerance = 1e-12
      )
      expect_equal(
        cost_point(bh, Y)$pcost,
        sum(cost_point(dense, Y)$pcost),
        tolerance = 1e-12
      )
    })
  }
})

test_that("exact Barnes-Hut Jensen-Shannon SNE matches dense on the saved graph", {
  fixture <- saved_neighbors_fixture()
  Y <- cbind(sin(seq_len(24)), cos(seq_len(24))) / 3
  Y[2, ] <- Y[1, ]
  for (kernel in c("knn", "nngaussian")) {
    local({
      prepare <- function(X, bh, ...) {
        smallvis_run(
          X,
          method = "jssne",
          bh = bh,
          inp_kernel = kernel,
          perplexity = 3,
          theta = 0,
          eps = 1e-5,
          n_threads = 2,
          prepare_only = TRUE,
          epoch_callback = FALSE,
          verbose = FALSE,
          ...
        )$cost_fn
      }
      dense <- prepare(fixture$D, FALSE)
      local_mocked_bindings(get_nn = function(...) stop("neighbor search ran"))
      bh <- prepare(fixture$X, TRUE, neighbors = fixture$graph)
      expect_s4_class(bh$P, "dgCMatrix")
      expect_equal(as.matrix(bh$P), dense$P, tolerance = 1e-12)
      expect_equal(cost_grad(bh, Y)$G, cost_grad(dense, Y)$G, tolerance = 1e-12)
      expect_equal(
        cost_point(bh, Y)$pcost,
        sum(cost_point(dense, Y)$pcost),
        tolerance = 1e-12
      )
    })
  }
})

test_that("invalid saved graphs and incompatible options fail before search", {
  fixture <- saved_neighbors_fixture()
  local_mocked_bindings(get_nn = function(...) stop("neighbor search ran"))
  fit <- function(graph = fixture$graph, ...) {
    args <- utils::modifyList(
      list(
        X = fixture$X,
        bh = TRUE,
        perplexity = 3,
        neighbors = graph,
        epoch_callback = FALSE,
        verbose = FALSE
      ),
      list(...)
    )
    do.call(smallvis, args)
  }
  expect_error(
    fit(method = "hlsne", bh = FALSE),
    "selected kernel is 'gaussian'"
  )
  expect_error(
    fit(method = "mmds", bh = FALSE),
    "does not use precomputed neighbors"
  )
  expect_error(fit(nn = "exact"), "Leave nn unset")
  expect_error(fit(method = list("tee", nn = "exact")), "Leave nn unset")
  expect_error(fit(inp_kernel = "gaussian"), "neighbors requires inp_kernel")
  expect_error(fit(perplexity = c(2, 3)), "single numeric perplexity")
  expect_error(fit(perplexity = "idp"), "single numeric perplexity")
  expect_error(fit(perplexity = 4), "Perplexity 4 .* requires 12")
  expect_error(
    fit(inp_kernel = "knn", perplexity = 2.5),
    "whole-number perplexities"
  )

  graph <- fixture$graph
  graph$dist <- graph$dist[-1, ]
  expect_error(fit(graph), "matching dimensions")
  graph <- fixture$graph
  graph$idx[1, 2] <- 2.5
  expect_error(fit(graph), "one-based integer indices")
  graph$idx[1, 2] <- graph$idx[1, 3]
  expect_error(fit(graph), "duplicate indices in row 1")
  graph <- fixture$graph
  graph$dist[1, 2] <- Inf
  expect_error(fit(graph), "finite nonnegative distances")
  graph$dist[1, 2] <- graph$dist[1, 10]
  expect_error(fit(graph), "not sorted by distance in row 1")
})

test_that("saved neighbors cover every active stage before optimization starts", {
  fixture <- saved_neighbors_fixture()
  args <- list(
    X = fixture$X,
    method = list("tee", inp_kernel = "nngaussian", perplexity = 3),
    inp_kernel = "knn",
    bh = TRUE,
    neighbors = fixture$graph,
    max_iter = 5,
    epoch_callback = FALSE,
    verbose = FALSE
  )
  local({
    local_mocked_bindings(find_beta_knn_cpp = function(...) {
      stop("calibration ran")
    })
    expect_error(
      do.call(smallvis_perpstep, c(args, list(step_iter = 2))),
      "Perplexity 8 .* requires 23"
    )
    args$method[[1L]] <- "hlsne"
    args$bh <- FALSE
    expect_error(
      do.call(smallvis_perpstep, c(args, list(step_iter = 2))),
      "Perplexity 8 .* requires 23"
    )
  })

  local_mocked_bindings(get_nn = function(...) stop("neighbor search ran"))
  direct <- do.call(smallvis_perpstep, c(args, list(step_iter = 0)))
  expect_true(all(is.finite(direct)))
  args$method$inp_kernel <- "knn"
  stepped <- do.call(smallvis_perpstep, c(args, list(step_iter = 2)))
  expect_equal(dim(stepped), c(24L, 2L))
  expect_true(all(is.finite(stepped)))
})

test_that("saved smooth kNN neighborhoods cover the whole schedule", {
  fixture <- saved_neighbors_fixture()
  args <- list(
    X = fixture$X,
    method = list("tumap", inp_kernel = "skd", perplexity = 4.5),
    inp_kernel = "knn",
    bh = TRUE,
    neighbors = list(
      idx = fixture$graph$idx[, 1:4],
      dist = fixture$graph$dist[, 1:4]
    ),
    max_iter = 3,
    epoch_callback = FALSE,
    verbose = FALSE
  )
  local({
    local_mocked_bindings(smooth_knn_distances = function(...) {
      stop("calibration ran")
    })
    expect_error(
      do.call(smallvis_perpstep, c(args, list(step_iter = 1))),
      "Perplexity 8 .* requires 7"
    )
    args$method <- list("umap", perplexity = 4.5)
    args$inp_kernel <- NULL
    args$bh <- FALSE
    expect_error(
      do.call(smallvis_perpstep, c(args, list(step_iter = 1))),
      "Perplexity 8 .* requires 7"
    )
  })

  local_mocked_bindings(get_nn = function(...) stop("neighbor search ran"))
  direct <- do.call(smallvis_perpstep, c(args, list(step_iter = 0)))
  expect_true(all(is.finite(direct)))
  args$neighbors <- fixture$graph
  stepped <- do.call(smallvis_perpstep, c(args, list(step_iter = 1)))
  expect_true(all(is.finite(stepped)))
})
