test_that("string and list methods share controls and list options take precedence", {
  X <- matrix(sin(seq_len(144)^2), 24, 6)
  run <- function(method, ...) {
    smallvis(
      X,
      method = method,
      scale = FALSE,
      perplexity = 4,
      max_iter = 0,
      Y_init = "rand",
      ret_extra = "P",
      epoch_callback = FALSE,
      verbose = FALSE,
      ...
    )
  }
  gaussian <- run("tsne")
  knn <- run("tsne", inp_kernel = "knn")
  expect_false(isTRUE(all.equal(gaussian$P, knn$P)))
  expect_equal(knn$P, run(list("TSNE"), inp_kernel = "knn")$P)
  expect_equal(
    knn$P,
    run(list("tSnE", inp_kernel = "knn"), inp_kernel = "gaussian")$P
  )
  expect_equal(run("umap")$P, run(list("UMAP"))$P)
  expect_equal(run("asne")$P, run(list("ASNE"))$P)
  expect_equal(
    run("tsne", inp_kernel = "knn")$P,
    run(list("tsne", perplexity = 4), inp_kernel = "knn")$P
  )
  expect_equal(
    run(list("tsne", perplexity = 6))$P,
    smallvis(
      X,
      scale = FALSE,
      perplexity = 6,
      max_iter = 0,
      ret_extra = "P",
      epoch_callback = FALSE,
      verbose = FALSE
    )$P
  )
})

test_that("resolved backend controls reach native and neighbor routines", {
  calls <- 0L
  native <- tsne_grad_cpp
  local_mocked_bindings(tsne_grad_cpp = function(...) {
    calls <<- calls + 1L
    native(...)
  })
  X <- matrix(sin(seq_len(144)^2), 24, 6)
  for (method in list("tsne", list("TSNE"))) {
    smallvis(
      X,
      method = method,
      use_cpp = TRUE,
      perplexity = 4,
      max_iter = 1,
      epoch_callback = FALSE,
      verbose = FALSE
    )
  }
  expect_equal(calls, 2L)
  neighbors <- get_nn
  # Capture the selected search policy without doing approximate search.
  seen <- logical()
  local_mocked_bindings(get_nn = function(X, k, exact, n_threads, verbose) {
    seen <<- c(seen, exact)
    neighbors(X, k, exact = TRUE, n_threads, verbose)
  })
  smallvis(
    X,
    inp_kernel = "knn",
    nn = "approximate",
    perplexity = 4,
    max_iter = 0,
    epoch_callback = FALSE,
    verbose = FALSE
  )
  smallvis(
    X,
    method = list("tsne", nn = "exact"),
    nn = "approximate",
    inp_kernel = "knn",
    perplexity = 4,
    max_iter = 0,
    epoch_callback = FALSE,
    verbose = FALSE
  )
  expect_equal(seen, c(FALSE, TRUE))
})

test_that("both method forms select sparse Barnes-Hut and accept dense Gaussian input", {
  X <- matrix(sin(seq_len(144)^2), 24, 6)
  Y <- cbind(sin(1:24), cos(1:24))
  run <- function(method, bh, kernel) {
    smallvis(
      X,
      method = method,
      bh = bh,
      inp_kernel = kernel,
      nn = "exact",
      use_cpp = TRUE,
      scale = FALSE,
      perplexity = 4,
      Y_init = Y,
      theta = 0,
      max_iter = 0,
      ret_extra = "P",
      epoch_callback = FALSE,
      verbose = FALSE
    )
  }
  for (kernel in c("gaussian", "nngaussian", "knn")) {
    bh <- run("tsne", TRUE, kernel)
    expect_s4_class(bh$P, "dgCMatrix")
    expect_equal(bh$P, run(list("tsne"), TRUE, kernel)$P)
    if (kernel == "gaussian") {
      dense <- run("tsne", FALSE, kernel)
      expect_equal(as.matrix(bh$P), dense$P, tolerance = 1e-12)
      expect_equal(bh$costs, sum(dense$costs), tolerance = 1e-12)
    }
  }
})

test_that("custom constructors retain their defaults and receive explicit controls", {
  seen <- NULL
  custom <- function(perplexity, eps = 0.123) {
    seen <<- c(perplexity = perplexity, eps = eps)
    tsne(perplexity = perplexity, eps = eps)
  }
  for (explicit in c(FALSE, TRUE)) {
    args <- list(
      X = iris10,
      method = custom,
      perplexity = 4,
      max_iter = 0,
      epoch_callback = FALSE,
      verbose = FALSE
    )
    if (explicit) {
      args$eps <- 0.01
    }
    do.call(smallvis, args)
    expect_equal(unname(seen), c(4, if (explicit) 0.01 else 0.123))
  }
})

test_that("bad dispatch options fail before neighbor search", {
  local_mocked_bindings(get_nn = function(...) stop("neighbor search ran"))
  for (method in list(
    list("tsne", bogus = TRUE),
    stats::setNames(list("tsne", 1, 2), c("", "eps", "eps")),
    list("tsne", 1),
    list()
  )) {
    expect_error(
      smallvis(
        iris10,
        method = method,
        bh = TRUE,
        epoch_callback = FALSE,
        verbose = FALSE
      ),
      "method|Method"
    )
  }
  expect_error(
    smallvis(iris10, method = "umap", bh = TRUE),
    "Barnes-Hut is only supported"
  )
  expect_error(
    smallvis(iris10, method = "tumap", bh = TRUE, inp_kernel = "skd"),
    "nngaussian"
  )
  expect_error(
    smallvis(iris10, method = "tumap", bh = TRUE, eps = 0),
    "requires a finite positive eps"
  )
})

test_that("reported perplexity follows the resolved method options", {
  res <- smallvis(
    iris10,
    method = list("tsne", perplexity = 6),
    perplexity = 4,
    max_iter = 1,
    ret_extra = TRUE,
    epoch_callback = FALSE,
    verbose = FALSE
  )
  expect_equal(res$perplexity, 6)
})

test_that("optional distances use the resolved method controls", {
  seen <- NULL
  distances <- calc_d
  local_mocked_bindings(calc_d = function(X, use_cpp, n_threads) {
    seen <<- list(use_cpp = use_cpp, n_threads = n_threads)
    distances(X, use_cpp = use_cpp, n_threads = n_threads)
  })
  smallvis(
    iris10,
    method = list("tsne", use_cpp = TRUE, n_threads = 2),
    perplexity = 4,
    max_iter = 0,
    ret_extra = "DX",
    epoch_callback = FALSE,
    verbose = FALSE
  )
  expect_equal(seen, list(use_cpp = TRUE, n_threads = 2))
})

test_that("inherited and specialized methods honor backend and thread controls", {
  native <- find_beta_cpp
  distances <- calc_d2
  input_threads <- output_threads <- integer()
  output_backends <- logical()
  local_mocked_bindings(
    find_beta_cpp = function(..., n_threads = 1L) {
      input_threads <<- c(input_threads, n_threads)
      native(..., n_threads = n_threads)
    },
    calc_d2 = function(X, use_cpp = FALSE, n_threads = 1L) {
      output_backends <<- c(output_backends, use_cpp)
      output_threads <<- c(output_threads, n_threads)
      distances(X, use_cpp = use_cpp, n_threads = n_threads)
    }
  )
  cases <- list(
    list("hssne"),
    list("dhssne"),
    list("snerv"),
    list("absne", alpha = 0.75, lambda = 0.5),
    list("absne", alpha = 1, lambda = 0),
    list("absne", alpha = 2, lambda = 2),
    list("absne", alpha = 0, lambda = 1),
    list("absne", alpha = 0, lambda = 0),
    list("jssne"),
    list("chsne"),
    list("hlsne")
  )
  # Equal input distances give the same normalized P in both calibrators even
  # though their precision initializations differ. Backend calls still matter.
  X <- diag(8)
  Y <- cbind(sin(1:8), cos(1:8))
  for (method in cases) {
    reference <- NULL
    for (nested in c(FALSE, TRUE)) {
      threads <- if (nested) 2L else 1L
      for (use_cpp in c(FALSE, TRUE)) {
        input_threads <- output_threads <- integer()
        output_backends <- logical()
        selector <- if (nested) {
          c(method, list(use_cpp = use_cpp, n_threads = threads))
        } else {
          method
        }
        fit <- smallvis(
          X,
          method = selector,
          use_cpp = if (nested) !use_cpp else use_cpp,
          n_threads = if (nested) 3L else threads,
          scale = FALSE,
          perplexity = 7,
          Y_init = Y,
          max_iter = 2,
          opt = list("adam", eta = 0.001),
          epoch = 1,
          tol = 0,
          min_cost = -Inf,
          ret_extra = "P",
          epoch_callback = FALSE,
          verbose = FALSE
        )
        expect_equal(input_threads, if (use_cpp) threads else integer())
        expect_gt(length(output_threads), 0)
        expect_equal(output_backends, rep(use_cpp, length(output_backends)))
        expect_equal(output_threads, rep(threads, length(output_threads)))
        if (is.null(reference)) {
          reference <- fit
        } else {
          expect_equal(fit$P, reference$P, tolerance = 1e-12)
          expect_equal(fit$Y, reference$Y, tolerance = 1e-10)
          expect_equal(fit$costs, reference$costs, tolerance = 1e-10)
        }
      }
    }
  }
})
