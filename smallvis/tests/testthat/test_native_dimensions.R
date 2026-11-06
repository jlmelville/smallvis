test_that("dense accelerated routes preserve general output dimensions", {
  X <- as.matrix(iris[1:12, 1:4])
  for (method in list("mmds", list("tsne", inp_kernel = "sigma"))) {
    for (k in c(1L, 2L, 3L)) {
      Y <- matrix(sin(seq_len(12 * k)), 12, k)
      run <- function(use_cpp) {
        current_method <- method
        if (is.list(current_method)) {
          current_method$use_cpp <- use_cpp
          current_method$n_threads <- 2
        }
        smallvis(
          X,
          k = k,
          method = method,
          Y_init = Y,
          perplexity = 1,
          max_iter = 1,
          eta = 0.01,
          epoch_callback = FALSE,
          verbose = FALSE,
          use_cpp = use_cpp,
          n_threads = 2
        )
      }
      expect_equal(dim(run(TRUE)), c(12L, k))
      expect_equal(run(TRUE), run(FALSE), tolerance = 1e-12)
    }
  }
})

test_that("Barnes-Hut rejects unsupported dimensions before calibration", {
  for (k in c(1, 3)) {
    expect_error(
      smallvis(
        iris10,
        bh = TRUE,
        k = k,
        perplexity = 30,
        epoch_callback = FALSE,
        verbose = FALSE
      ),
      "Barnes-Hut requires k = 2"
    )
  }
  # Guard each native entry because coordinates can also come from optimizers.
  for (k in c(1, 3)) {
    Y <- matrix(0, 3, k)
    expect_error(
      bh_plogq_cpp(integer(), rep(0L, 4), numeric(), Y),
      "two output dimensions"
    )
    expect_error(
      bh_tsne_gradient_cpp(integer(), rep(0L, 4), numeric(), Y),
      "two output dimensions"
    )
    expect_error(tweight_cpp(Y), "two output dimensions")
    expect_error(
      mmds_grad_cpp(matrix(0, 3, 3), matrix(0, 3, 3), Y, 1e-16, 1L),
      "two output dimensions"
    )
    expect_error(
      tsne_grad_cpp(matrix(0, 3, 3), matrix(0, 3, 3), 1, Y, 1L),
      "two output dimensions"
    )
  }
})
