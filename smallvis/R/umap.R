# UMAP --------------------------------------------------------------------

umap <- function(
  perplexity,
  inp_kernel = "skd",
  symmetrize = "umap",
  spread = 1,
  min_dist = 0.001,
  gr_eps = 0.1,
  eps = 1e-9,
  row_weight = NULL,
  n_threads = 0,
  use_cpp = FALSE
) {
  if (!is.null(row_weight)) {
    row_normalize <- row_weight
  } else {
    row_normalize <- FALSE
  }
  lreplace(
    tsne(perplexity, use_cpp = use_cpp, n_threads = n_threads),
    gradient_only = gr_eps != 0,
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      symmetrize <- match.arg(tolower(symmetrize), true_symmetrize_options())
      cost <- sne_init(
        cost,
        X,
        perplexity = perplexity,
        kernel = inp_kernel,
        symmetrize = symmetrize,
        normalize = FALSE,
        row_normalize = row_normalize,
        n_threads = n_threads,
        verbose = verbose,
        ret_extra = ret_extra,
        use_cpp = use_cpp
      )

      cost <- init_ab(
        cost,
        spread = spread,
        min_dist = min_dist,
        verbose = verbose
      )

      cost$eps <- eps
      cost
    },
    cache_input = function(cost) {
      P <- cost$P
      eps <- cost$eps
      cost$Cp <- colSums(P * log(P + eps) + (1 - P) * log1p(-P + eps))
      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)

      P <- cost$P
      eps <- cost$eps
      W <- cost$W
      cost$pcost <- colSums(-P * log(W + eps) - (1 - P) * log1p(-W + eps)) +
        cost$Cp
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)

      cost$G <- k2g(
        Y,
        4 * (cost$b / (cost$D2 + cost$eps + gr_eps)) * (cost$P - cost$W)
      )
      cost
    },
    update = function(cost, Y) {
      D2 <- calc_d2(Y, use_cpp = use_cpp, n_threads = n_threads)
      D2[D2 < 0] <- 0

      W <- 1 / (1 + cost$a * D2^cost$b)
      diag(W) <- 0

      cost$W <- W
      cost$D2 <- D2
      cost
    },
    export = cost_export
  )
}

# UMAP with the output kernel fixed to the t-distribution
tumap <- function(
  perplexity,
  inp_kernel = "skd",
  symmetrize = "umap",
  gr_eps = 0.1,
  eps = 1e-9,
  row_weight = NULL,
  n_threads = 0,
  use_cpp = FALSE
) {
  if (!is.null(row_weight)) {
    row_normalize <- row_weight
  } else {
    row_normalize <- FALSE
  }
  lreplace(
    umap(perplexity, n_threads = n_threads, use_cpp = use_cpp),
    gradient_only = gr_eps != 0,
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      symmetrize <- match.arg(tolower(symmetrize), true_symmetrize_options())
      cost$eps <- eps

      cost <- sne_init(
        cost,
        X,
        perplexity = perplexity,
        kernel = inp_kernel,
        symmetrize = symmetrize,
        normalize = FALSE,
        row_normalize = row_normalize,
        n_threads = n_threads,
        verbose = verbose,
        ret_extra = ret_extra,
        use_cpp = use_cpp
      )
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)

      cost$G <- k2g(
        Y,
        4 * (cost$W / ((1 - cost$W) + cost$eps + gr_eps)) * (cost$P - cost$W)
      )
      cost
    },
    update = function(cost, Y) {
      W <- calc_d2(Y, use_cpp = use_cpp, n_threads = n_threads)
      W <- 1 / (1 + W)
      diag(W) <- 0

      cost$W <- W
      cost
    }
  )
}

# t-UMAP where output and input affinities are normalized
ntumap <- function(
  perplexity,
  inp_kernel = "skd",
  symmetrize = "umap",
  gr_eps = 0.1,
  eps = 1e-9,
  n_threads = 0,
  use_cpp = FALSE
) {
  lreplace(
    tumap(perplexity, use_cpp = use_cpp, n_threads = n_threads),
    gradient_only = FALSE,
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      symmetrize <- match.arg(tolower(symmetrize), true_symmetrize_options())
      cost$eps <- eps
      cost <- sne_init(
        cost,
        X,
        perplexity = perplexity,
        kernel = inp_kernel,
        symmetrize = symmetrize,
        normalize = TRUE,
        n_threads = n_threads,
        verbose = verbose,
        ret_extra = ret_extra,
        use_cpp = use_cpp
      )
      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)

      P <- cost$P
      eps <- cost$eps
      Q <- cost$Q

      cost$pcost <- colSums(-P * logm(Q, eps) - (1 - P) * log1p(-Q + eps)) +
        cost$Cp
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      cost$G <- k2g(Y, 4 * cost$W * (cost$C - cost$sumC * cost$Q))
      cost
    },
    update = function(cost, Y) {
      P <- cost$P
      W <- calc_d2(Y, use_cpp = use_cpp, n_threads = n_threads)
      W <- 1 / (1 + W)
      diag(W) <- 0

      Q <- W / sum(W)
      C <- (P - Q) / (1 - Q)
      sumC <- sum(C)

      cost$W <- W
      cost$Q <- Q
      cost$C <- C
      cost$sumC <- sumC

      cost
    },
    export = cost_export
  )
}

# Fits a kernel for the output distances of the form w = 1 / (1 + a dsq ^ b)
# where dsq is the squared Euclidean distance.
# Standard t-SNE function is a = 1, b = 1.
# Default UMAP values are a = 1.929, b = 0.7915.
find_ab_params <- function(spread = 1, min_dist = 0.001) {
  xv <- seq(
    from = 0,
    to = spread * 3,
    length.out = 300
  )
  yv <- rep(0, length(xv))
  yv[xv < min_dist] <- 1
  yv[xv >= min_dist] <- exp(-(xv[xv >= min_dist] - min_dist) / spread)
  result <- try(
    {
      stats::nls(
        yv ~ 1 / (1 + a * xv^(2 * b)),
        start = list(a = 1, b = 1)
      )$m$getPars()
    },
    silent = TRUE
  )
  if (methods::is(result, "try-error")) {
    stop("Can't find a, b for provided spread/min_dist values")
  }
  result
}

# The UMAP equivalent of perplexity calibration in x2aff. k is continuous rather
# than integral and so is analogous to perplexity.
# Some differences:
# 1. The target value is the log2 of k, not the Shannon entropy associated
# with the desired perplexity.
# 2. Input weights are exponential, rather than Gaussian, with respect to the
# distances. The distances are also centered with respect to the smoothed
# distance to the nearest (non-zero distance) neighbor. A non-integral
# 'local_connectivity' value can result in this shortest distance between an
# interpolated value between two distances.
# 3. The weights are not normalized. Their raw sum is compared to the target
# value.
# 4. Distances beyond the k-nearest neighbors are not used in the calibration.
# The equivalent weights are set to 0.
# 5. Weights associated with distances shorter than the smoothed nearest
# neighbor distance are clamped to 1.
# This code has been converted from the original Python and may not be very
# idiomatic (or vectorizable).
# tol is SMOOTH_K_TOLERANCE in the Python code.
smooth_knn_distances <-
  function(
    X,
    k,
    n_iter = 64,
    local_connectivity = 1.0,
    bandwidth = 1.0,
    tol = 1e-5,
    min_k_dist_scale = 1e-3,
    cardinality = log2(k),
    n_threads = 0,
    verbose = FALSE,
    neighbors = NULL,
    sparse = FALSE,
    exact = FALSE
  ) {
    tsmessage("Commencing smooth kNN distance calibration for k = ", formatC(k))

    n <- if (methods::is(X, "dist")) attr(X, "Size") else nrow(X)
    validate_skd_k(k, n)
    # k controls a continuous target; only the neighbor count is rounded down.
    n_neighbors <- floor(k) - 1L
    if (is.null(neighbors)) {
      neighbors <- get_nn(
        X,
        k = n_neighbors + 1L,
        exact = exact,
        n_threads = n_threads,
        verbose = verbose
      )
    }
    nn_idx <- matrix(0L, nrow = n, ncol = n_neighbors)
    nn_dist <- matrix(0, nrow = n, ncol = n_neighbors + 1L)
    for (i in seq_len(n)) {
      keep <- utils::head(which(neighbors$idx[i, ] != i), n_neighbors)
      nn_idx[i, ] <- neighbors$idx[i, keep]
      nn_dist[i, -1] <- neighbors$dist[i, keep]
    }

    target <- cardinality * bandwidth
    rho <- rep(0, n)
    sigma <- rep(0, n)
    P <- matrix(0, nrow = n, ncol = if (sparse) n_neighbors else n)
    mean_distances <- NULL

    for (i in 1:n) {
      lo <- 0.0
      hi <- Inf
      mid <- 1.0

      ith_distances <- nn_dist[i, ]
      non_zero_dists <- ith_distances[ith_distances > 0.0]
      if (length(non_zero_dists) >= local_connectivity) {
        index <- floor(local_connectivity)
        interpolation <- local_connectivity - index
        if (index > 0) {
          if (interpolation <= tol) {
            rho[i] <- non_zero_dists[index]
          } else {
            rho[i] <- non_zero_dists[index] +
              interpolation *
                (non_zero_dists[index + 1] - non_zero_dists[index])
          }
        } else {
          rho[i] <- interpolation * non_zero_dists[1]
        }
      } else if (length(non_zero_dists) > 0) {
        rho[i] <- max(non_zero_dists)
      } else {
        rho[i] <- 0.0
      }

      for (iter in 1:n_iter) {
        psum <- 0.0
        for (j in 2:ncol(nn_dist)) {
          dist <- max(0, (nn_dist[i, j] - rho[i]))
          psum <- psum + exp(-(dist / mid))
        }
        val <- psum

        if (abs(val - target) < tol) {
          break
        }

        if (val > target) {
          hi <- mid
          mid <- (lo + hi) / 2.0
        } else {
          lo <- mid
          if (is.infinite(hi)) {
            mid <- mid * 2
          } else {
            mid <- (lo + hi) / 2.0
          }
        }
      }
      sigma[i] <- mid

      if (rho[i] > 0.0) {
        sigma[i] <- max(sigma[i], min_k_dist_scale * mean(ith_distances))
      } else {
        if (is.null(mean_distances)) {
          mean_distances <- mean(nn_dist)
        }
        sigma[i] <- max(sigma[i], min_k_dist_scale * mean_distances)
      }

      prow <- exp(-(nn_dist[i, ] - rho[i]) / (sigma[i] * bandwidth))
      prow[nn_dist[i, ] - rho[i] <= 0] <- 1
      P[i, if (sparse) seq_len(n_neighbors) else nn_idx[i, ]] <- prow[-1L]
    }
    if (sparse) {
      P <- Matrix::sparseMatrix(
        i = rep(seq_len(n), times = n_neighbors),
        j = as.vector(nn_idx),
        x = as.vector(P),
        dims = c(n, n)
      )
    }

    if (verbose) {
      summarize(sigma, "sigma summary", verbose = verbose)
    }
    list(sigma = sigma, rho = rho, P = P)
  }

validate_skd_k <- function(k, n) {
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k < 2 || k > n) {
    stop("UMAP k must be finite and between 2 and the number of observations")
  }
}


# set_op_mix_ratio = between 0 and 1 mixes in fuzzy set intersection
# set to 0 for intersection only
fuzzy_set_union <- function(X, set_op_mix_ratio = 1) {
  Xt <- Matrix::t(X)
  XX <- X * Xt
  set_op_mix_ratio * (X + Xt - XX) + (1 - set_op_mix_ratio) * XX
}

init_ab <- function(cost, spread = 1, min_dist = 0.001, verbose = FALSE) {
  ab_params <- find_ab_params(spread = spread, min_dist = min_dist)
  a <- ab_params[1]
  b <- ab_params[2]
  if (verbose) {
    message("Umap curve parameters = ", formatC(a), ", ", formatC(b))
  }
  cost$a <- a
  cost$b <- b
  cost
}
