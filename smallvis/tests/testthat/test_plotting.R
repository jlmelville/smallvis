test_that("the default callback plots with the public vizier interface", {
  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit(
    {
      grDevices::dev.off()
      unlink(plot_file)
    },
    add = TRUE
  )

  for (X in list(iris[1:10, ], iris10, stats::dist(iris10))) {
    res <- smallvis(
      X,
      Y_init = iris10_Y,
      perplexity = 5,
      max_iter = 2,
      epoch = 1,
      n_threads = 1,
      verbose = FALSE
    )

    expect_equal(dim(res), c(10L, 2L), info = class(X)[1])
  }
})
