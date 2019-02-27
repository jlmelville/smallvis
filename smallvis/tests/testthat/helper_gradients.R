gradient_fd <- function(Y, cost_fn, diff = .Machine$double.eps ^ (1 / 3)) {
  cost <- function(Y) {
    cost_fn <- cost_clear(cost_fn)
    cost_fn <- cost_fn$gr(cost_fn, Y)
    cost_fn <- cost_fn$pfn(cost_fn, Y)
    sum(cost_fn$pcost)
  }

  nr <- nrow(Y)
  nc <- ncol(Y)

  Gfd <- matrix(0, nrow = nr, ncol = nc)
  for (i in 1:nr) {
    for (j in 1:nc) {
      Yij_old <- Y[i, j]

      Y[i, j] <- Yij_old + diff
      cost_fwd <- cost(Y)

      Y[i, j] <- Yij_old - diff
      cost_back <- cost(Y)

      fd <- (cost_fwd - cost_back) / (2 * diff)
      Gfd[i, j] <- fd

      Y[i, j] <- Yij_old
    }
  }

  Gfd
}

expect_grad <- function(cost_fn,
                        Y,
                        label = "",
                        info = label,
                        diff = .Machine$double.eps ^ (1 / 3),
                        tolerance = 1e-6,
                        scale = 1) {

  cost_fn <- cost_clear(cost_fn)
  cost_fn <- cost_grad(cost_fn, Y)
  gan <- cost_fn$G
  gfd <- gradient_fd(Y, cost_fn)
  
  expect_false(any(is.nan(gan)))
  expect_false(any(is.nan(gfd)))
  expect_equal(gan, gfd, tolerance = tolerance, scale = scale,
               label = label, info = info,
               expected.label = "finite difference gradient")
}

test_grad <- function(method_name, tolerance = 1e-6, X = iris10, Y = iris10_Y,
                      ...) {
  if (!exists(method_name)) {
    return(testthat::fail(paste0("No method called '", method_name, "'")))
  }

  args <- list(...)
  cost_fn <- do.call(get(method_name), args)
  cost_fn <- cost_init(cost_fn, X, verbose = FALSE)
  expect_grad(cost_fn, Y, tolerance = tolerance,
              label = paste0(method_name, " ", l2s(args)))
}

l2s <- function(l) {
  if (length(l) == 0) {
    return("")
  }
  paste0(mapply(FUN = paste0, names(l), ":", l), collapse = " ")
}

# extract the last evaluated error
final_cost <- function(res, digits = 4) {
  signif(as.numeric(res$itercosts[length(res$itercosts)]), digits = digits)
}

# Covert a vector into a 2D matrix for generating Y output
c2y <- function(...) {
  matrix(unlist(list(...)), ncol = 2)
}
