# Generic Weight/Cost/Gradient functions ----------------------------------

logm <- function(m, eps = .Machine$double.eps) {
  diag(m) <- eps
  m <- log(m)
  diag(m) <- 0
  m
}

divm <- function(m, n, eps = .Machine$double.eps) {
  diag(m) <- eps
  m <- m / n
  diag(m) <- 0
  m
}

powm <- function(m, n, eps = .Machine$double.eps) {
  diag(m) <- eps
  m <- m^n
  diag(m) <- 0
  m
}

# Convert Force constant to Gradient
k2g <- function(Y, K, symmetrize = FALSE) {
  if (symmetrize) {
    K <- K + t(K)
  }
  Y * colSums(K) - (K %*% Y)
}

cost_init <- function(cost,
                      X,
                      max_iter,
                      verbose = FALSE,
                      ret_extra = c()) {
  if (!is.null(cost$init)) {
    cost <- cost$init(
      cost,
      X,
      verbose = verbose,
      ret_extra = ret_extra,
      max_iter = max_iter
    )
  }
  cost_cache_input(cost)
}

cost_grad <- function(cost, Y) {
  cost$gr(cost, Y)
}

cost_point <- function(cost, Y) {
  cost$pfn(cost, Y)
}

# Clear values dependent on Y (or other parameters)
# This is the mechanism by which Gradient and Function evaluations will
# detect that they need recalculate distances, weights, probabilities etc.
cost_clear <- function(cost) {
  if (!is.null(cost$sentinel)) {
    cost[[cost$sentinel]] <- NULL
  } else if (is.null(cost$clear)) {
    cost <- cost$clear(cost)
  }
  cost
}

# Update matrices dependent on Y
# Check for a NULL value of the sentinel to avoid unnecessary recalculation
cost_update <- function(cost, Y) {
  if (!is.null(cost$sentinel)) {
    if (is.null(cost[[cost$sentinel]])) {
      cost <- cost$update(cost, Y)
    }
  } else {
    cost <- cost$update(cost, Y)
  }
  cost
}

# Calculate the cost function value. If opt_res$f is non-NULL, use that
cost_eval <- function(cost, Y, opt_res = NULL) {
  if (is.null(opt_res) || is.null(opt_res$f)) {
    cost <- cost_point(cost, Y)
    pcosts <- cost$pcost
    cost_val <- sum(pcosts)
  } else {
    cost_val <- opt_res$f
  }

  list(cost = cost, value = cost_val)
}

# Default export of values associated with a method
# If the value is associated with the cost list, e.g. cost$P, then asking
# for val = "P" or val = "p" will return it. Otherwise, returns NULL
cost_export <- function(cost, val) {
  res <- NULL
  if (!is.null(cost[[val]])) {
    res <- cost[[val]]
  } else if (!is.null(cost[[toupper(val)]])) {
    res <- cost[[toupper(val)]]
  }
  res
}

cost_cache_input <- function(cost) {
  if (!is.null(cost$cache_input)) {
    cost <- cost$cache_input(cost)
  }
  cost
}

start_exaggerating <- function(cost, exaggeration_factor) {
  if (!is.null(cost$exaggerate) && exaggeration_factor != 1) {
    cost <- cost$exaggerate(cost, exaggeration_factor)
  }
  cost_cache_input(cost)
}

stop_exaggerating <- function(cost, exaggeration_factor) {
  if (!is.null(cost$exaggerate) && exaggeration_factor != 1) {
    cost <- cost$exaggerate(cost, 1 / exaggeration_factor)
  }
  cost_cache_input(cost)
}
