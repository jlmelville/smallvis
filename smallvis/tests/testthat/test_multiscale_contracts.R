test_that("multiscale consumers derive their inputs from the averaged P", {
  X <- matrix(sin(seq_len(144)^2), 24, 6)
  Y <- cbind(sin(1:24), cos(1:24))
  for (method in c("tsne", "asne", "wtsne", "wssne", "tsneu", "pstsne")) {
    run <- function(perplexity, iterations = 2) {
      smallvis(
        X,
        method = method,
        perplexity = perplexity,
        scale = FALSE,
        Y_init = Y,
        n_threads = 1,
        max_iter = iterations,
        opt = list("adam", eta = 0.001),
        epoch = 1,
        tol = 0,
        min_cost = -Inf,
        ret_extra = c("P", "pdeg", "M", "V"),
        epoch_callback = FALSE,
        verbose = FALSE
      )
    }
    single <- run(4)
    singleton <- run(list("multiscale", 4))
    mixture <- run(list("multiscale", c(4, 6)))
    expect_equal(singleton$P, single$P, tolerance = 1e-10)
    expect_equal(singleton$Y, single$Y, tolerance = 1e-10)
    expect_equal(singleton$costs, single$costs, tolerance = 1e-10)
    expect_equal(singleton$itercosts, single$itercosts, tolerance = 1e-10)
    expect_equal(
      mixture$P,
      (single$P + run(6, 0)$P) / 2,
      tolerance = 1e-10
    )
    for (fit in list(singleton, mixture)) {
      expect_equal(fit$iter, 2)
      expect_true(all(is.finite(fit$Y)))
      expect_true(all(is.finite(fit$costs)))
      expect_equal(diag(fit$P), rep(0, 24))
      if (method == "asne") {
        expect_equal(rowSums(fit$P), rep(1, 24))
      } else {
        expect_equal(sum(fit$P), 1)
      }
      expect_equal(fit$pdeg, colSums(fit$P))
      if (method %in% c("wtsne", "wssne")) {
        degrees <- colSums(fit$P)
        expect_equal(fit$M, outer(degrees, degrees))
      }
      if (method %in% c("tsneu", "pstsne")) {
        expect_equal(fit$V, 24 * fit$P, tolerance = 1e-10)
      }
    }
  }
})

test_that("automatic bandwidth transfer rejects multiscale before calibration", {
  local_mocked_bindings(msp = function(...) stop("calibration ran"))
  X <- matrix(sin(seq_len(96)), 24, 4)
  for (method in c("btsne", "bssne", "basne", "btasne", "bnerv")) {
    for (nested in c(FALSE, TRUE)) {
      selector <- if (nested) {
        list(method, perplexity = list("multiscale", c(4, 6)))
      } else {
        method
      }
      expect_error(
        smallvis(
          X,
          method = selector,
          perplexity = if (nested) 4 else list("multiscale", c(4, 6)),
          max_iter = 1,
          epoch_callback = FALSE,
          verbose = FALSE
        ),
        paste0(method, ".*multiscale")
      )
    }
  }
})

test_that("explicit output beta remains compatible with multiscale affinities", {
  X <- matrix(sin(seq_len(144)^2), 24, 6)
  Y <- cbind(sin(1:24), cos(1:24))
  for (method in c("btsne", "bssne", "basne", "btasne")) {
    run <- function(perplexity) {
      smallvis(
        X,
        method = list(method, beta = 0.5),
        perplexity = perplexity,
        scale = FALSE,
        Y_init = Y,
        max_iter = 2,
        opt = list("adam", eta = 0.001),
        epoch = 1,
        tol = 0,
        min_cost = -Inf,
        ret_extra = c("P", "beta"),
        epoch_callback = FALSE,
        verbose = FALSE
      )
    }
    single <- run(4)
    singleton <- run(list("multiscale", 4))
    mixture <- run(list("multiscale", c(4, 6)))
    expect_equal(singleton$P, single$P, tolerance = 1e-10)
    expect_equal(singleton$Y, single$Y, tolerance = 1e-10)
    expect_equal(singleton$costs, single$costs, tolerance = 1e-10)
    expect_equal(mixture$beta, 0.5)
    expect_equal(mixture$P, (single$P + run(6)$P) / 2, tolerance = 1e-10)
    expect_equal(mixture$iter, 2)
    expect_true(all(is.finite(mixture$Y)))
    expect_true(all(is.finite(mixture$costs)))
  }
})
