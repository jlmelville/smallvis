# UMAP --------------------------------------------------------------------

umap <- function(perplexity,
                 inp_kernel = "skd",
                 symmetrize = "umap",
                 spread = 1,
                 min_dist = 0.001,
                 gr_eps = 0.1,
                 eps = 1e-9,
                 row_weight = NULL,
                 n_threads = 0,
                 use_cpp = FALSE) {
  if (!is.null(row_weight)) {
    row_normalize <- row_weight
  } else {
    row_normalize <- FALSE
  }
  lreplace(
    tsne(perplexity, use_cpp = use_cpp, n_threads = n_threads),
    init = function(cost,
                    X,
                    max_iter,
                    verbose = FALSE,
                    ret_extra = c()) {
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

      cost <- init_ab(cost,
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
      cost$pcost <- colSums(-P * log(W + eps) - (1 - P) * log1p(-W + eps)) + cost$Cp
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)

      cost$G <- k2g(Y, 4 * (cost$b / (cost$D2 + cost$eps + gr_eps)) * (cost$P - cost$W))
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
tumap <- function(perplexity,
                  inp_kernel = "skd",
                  symmetrize = "umap",
                  gr_eps = 0.1,
                  eps = 1e-9,
                  row_weight = NULL,
                  n_threads = 0,
                  use_cpp = FALSE) {
  if (!is.null(row_weight)) {
    row_normalize <- row_weight
  } else {
    row_normalize <- FALSE
  }
  lreplace(
    umap(perplexity, n_threads = n_threads, use_cpp = use_cpp),
    init = function(cost,
                    X,
                    max_iter,
                    verbose = FALSE,
                    ret_extra = c()) {
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

      cost$G <- k2g(Y, 4 * (cost$W / ((1 - cost$W) + cost$eps + gr_eps)) * (cost$P - cost$W))
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
ntumap <- function(perplexity,
                   inp_kernel = "skd",
                   symmetrize = "umap",
                   gr_eps = 0.1,
                   eps = 1e-9,
                   n_threads = 0,
                   use_cpp = FALSE) {
  lreplace(
    tumap(perplexity, use_cpp = use_cpp, n_threads = n_threads),
    init = function(cost,
                    X,
                    max_iter,
                    verbose = FALSE,
                    ret_extra = c()) {
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

      cost$pcost <- colSums(-P * logm(Q, eps) - (1 - P) * log1p(-Q + eps)) + cost$Cp
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
      stats::nls(yv ~ 1 / (1 + a * xv^(2 * b)), start = list(a = 1, b = 1))$m$getPars()
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
  function(X,
           k,
           n_iter = 64,
           local_connectivity = 1.0,
           bandwidth = 1.0,
           tol = 1e-5,
           min_k_dist_scale = 1e-3,
           cardinality = log2(k),
           n_threads = 0,
           verbose = FALSE) {
    tsmessage("Commencing smooth kNN distance calibration for k = ", formatC(k))

    if (methods::is(X, "dist")) {
      X <- as.matrix(X)
      nn_idx <- t(apply(X, 2, order))[, 1:k]
      nn_dist <- matrix(0, nrow = nrow(X), ncol = k)
      for (i in 1:nrow(nn_idx)) {
        nn_dist[i, ] <- X[i, nn_idx[i, ]]
      }
    } else {
      # TODO: shouldn't this be k + 1 ?
      knn <- get_nn(X, k = k, n_threads = n_threads, verbose = verbose)
      knn$idx <- knn$idx[, 2:k]
      knn$dist <- knn$dist[, 2:k]

      nn_idx <- matrix(nrow = nrow(X), ncol = k)
      nn_idx[, 1] <- 1:nrow(nn_idx)
      nn_idx[, 2:ncol(nn_idx)] <- knn$idx
      nn_dist <- matrix(0, nrow = nrow(X), ncol = k)
      nn_dist[, 2:ncol(nn_dist)] <- knn$dist
    }

    n <- nrow(nn_dist)
    target <- cardinality * bandwidth
    rho <- rep(0, n)
    sigma <- rep(0, n)
    P <- matrix(0, nrow = n, ncol = n)
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
            rho[i] <- non_zero_dists[index] + interpolation *
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
      P[i, nn_idx[i, ]] <- prow
    }
    diag(P) <- 0

    if (verbose) {
      summarize(sigma, "sigma summary", verbose = verbose)
    }
    list(sigma = sigma, rho = rho, P = P)
  }


# set_op_mix_ratio = between 0 and 1 mixes in fuzzy set intersection
# set to 0 for intersection only
fuzzy_set_union <- function(X, set_op_mix_ratio = 1) {
  XX <- X * t(X)
  set_op_mix_ratio * (X + t(X) - XX) + (1 - set_op_mix_ratio) * XX
}

init_ab <- function(cost,
                    spread = 1,
                    min_dist = 0.001,
                    verbose = FALSE) {
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
