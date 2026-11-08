# Generic Weight/Cost/Gradient functions ----------------------------------

# Convert Force constant to Gradient
k2g <- function(Y, K, symmetrize = FALSE) {
  if (symmetrize) {
    K <- K + t(K)
  }
  Y * colSums(K) - (K %*% Y)
}

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

# Shift log weights before exponentiating. Conditional distributions use
# independent row shifts; a globally normalized matrix uses one shared shift.
exp_shift <- function(X, matrix_normalize = FALSE) {
  diag(X) <- -Inf
  if (matrix_normalize) {
    exp(X - max(X))
  } else {
    exp(sweep(X, 1, apply(X, 1, max)))
  }
}

expQ <- function(
  Y,
  eps = .Machine$double.eps,
  beta = NULL,
  A = NULL,
  is_symmetric = FALSE,
  matrix_normalize = FALSE,
  use_cpp = FALSE,
  n_threads = 1
) {
  logW <- -calc_d2(Y, use_cpp = use_cpp, n_threads = n_threads)
  if (!is.null(beta)) {
    logW <- logW * beta
  }
  if (!is.null(A)) {
    logW <- logW + log(A)
  }
  W <- exp_shift(logW, matrix_normalize = matrix_normalize)
  if (matrix_normalize) {
    Z <- sum(W)
  } else {
    # Row shifts break symmetry even when the input affinities are symmetric.
    Z <- rowSums(W)
  }
  # cost of division (vs storing 1/Z and multiplying) seems small
  Q <- W / Z

  if (eps > 0) {
    Q[Q < eps] <- eps
  }
  diag(Q) <- 0

  list(Q = Q, Z = Z)
}

# KL divergence using Q directly
kl_costQ <- function(cost, Y) {
  cost <- cost_update(cost, Y)

  # P log(P / Q) = P log P - P log Q
  cost$pcost <- cost$plogp - colSums(cost$P * logm(cost$Q, cost$eps))
  cost
}

kl_costQr <- function(cost, Y) {
  cost <- cost_update(cost, Y)

  # P log(P / Q) = P log P - P log Q
  cost$pcost <- cost$plogp - rowSums(cost$P * logm(cost$Q, cost$eps))
  cost
}

kl_cost <- function(cost, Y) {
  cost <- cost_update(cost, Y)
  # P log(P / Q) = P log P - P log Q
  cost$pcost <- cost$plogp - colSums(cost$P * logm(cost$W / cost$Z, cost$eps))
  cost
}

kl_costr <- function(cost, Y) {
  cost <- cost_update(cost, Y)
  cost$pcost <- cost$plogp - rowSums(cost$P * logm(cost$W / cost$Z, cost$eps))
  cost
}

cost_init <- function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
  has_neighbors <- !is.null(cost[["neighbors"]])
  if (!is.null(cost$init)) {
    cost <- cost$init(
      cost,
      X,
      verbose = verbose,
      ret_extra = ret_extra,
      max_iter = max_iter
    )
  }
  if (has_neighbors) {
    if (!isTRUE(cost[["neighbors_used"]])) {
      stop("The selected method does not use precomputed neighbors")
    }
    cost$neighbors_used <- NULL
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
  } else if (!is.null(cost$clear)) {
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
