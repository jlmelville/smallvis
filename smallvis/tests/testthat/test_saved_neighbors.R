saved_neighbors_fixture <- function() {
  X <- cbind(seq_len(24), sin(seq_len(24)), cos(seq_len(24)))
  D <- stats::dist(X)
  list(X = X, D = D, graph = get_nn(D, k = 10L))
}

test_that("all Barnes-Hut methods can calibrate the supplied graph without search", {
  fixture <- saved_neighbors_fixture()
  Y <- cbind(sin(seq_len(24)), cos(seq_len(24))) / 3
  run <- function(method, kernel, X, ...) {
    smallvis(
      X,
      method = method,
      bh = TRUE,
      inp_kernel = kernel,
      perplexity = 3,
      Y_init = Y,
      theta = 0,
      max_iter = 1,
      ret_extra = "P",
      epoch_callback = FALSE,
      verbose = FALSE,
      ...
    )
  }
  cases <- expand.grid(
    method = c("tsne", "hssne", "tee"),
    kernel = c("knn", "nngaussian"),
    stringsAsFactors = FALSE
  )
  expected <- lapply(seq_len(nrow(cases)), function(i) {
    run(cases$method[i], cases$kernel[i], fixture$D, nn = "exact")
  })
  local_mocked_bindings(get_nn = function(...) stop("neighbor search ran"))
  for (i in seq_len(nrow(cases))) {
    graph <- fixture$graph
    # Test both omitted self and self outside the first column.
    cols <- if (cases$kernel[i] == "knn") 2:10 else c(2:10, 1L)
    graph$idx <- graph$idx[, cols]
    graph$dist <- graph$dist[, cols]
    if (cases$method[i] == "tee") {
      storage.mode(graph$idx) <- "double"
    }
    actual <- run(
      cases$method[i],
      cases$kernel[i],
      fixture$X,
      neighbors = graph
    )
    expect_equal(actual$P, expected[[i]]$P)
    expect_equal(actual$Y, expected[[i]]$Y)
    expect_equal(actual$costs, expected[[i]]$costs)
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
  expect_error(fit(bh = FALSE), "neighbors requires Barnes-Hut")
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
    local_mocked_bindings(sne_init = function(...) stop("calibration ran"))
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
