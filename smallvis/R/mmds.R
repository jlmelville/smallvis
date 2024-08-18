# Distance Preserving Methods ---------------------------------------------

mmds_init <- function(cost,
                      X,
                      max_iter,
                      eps = .Machine$double.eps,
                      verbose = FALSE,
                      ret_extra = c(),
                      use_cpp = FALSE,
                      n_threads = 1) {
  tsmessage("Calculating pairwise distances")
  if (methods::is(X, "dist")) {
    cost$R <- as.matrix(X)
  } else {
    cost$R <- calc_d(X, use_cpp = use_cpp, n_threads = n_threads)
  }
  cost$eps <- eps
  cost
}

# Metric MDS, minimizing strain.
mmds <- function(eps = .Machine$double.eps,
                 use_cpp = FALSE,
                 n_threads = 1) {
  list(
    init = function(cost,
                    X,
                    max_iter,
                    verbose = FALSE,
                    ret_extra = c()) {
      cost <- mmds_init(
        cost,
        X,
        max_iter,
        eps,
        verbose,
        ret_extra,
        use_cpp = use_cpp,
        n_threads = n_threads
      )
      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)
      cost$pcost <- colSums((cost$R - cost$D)^2)
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      if (use_cpp) {
        cost$G <- mmds_grad_cpp(cost$R, cost$D, Y, eps = eps, n_threads = n_threads)
      } else {
        cost$G <- k2g(Y, -4 * (cost$R - cost$D) / (cost$D + cost$eps))
      }
      cost
    },
    update = function(cost, Y) {
      cost$D <- calc_d(Y, use_cpp = use_cpp, n_threads = n_threads)
      cost
    },
    sentinel = "D",
    export = function(cost, val) {
      res <- NULL
      switch(val,
        dx = {
          res <- cost$R
        },
        dy = {
          res <- cost$D
        }
      )
      res
    }
  )
}

smmds <- function(eps = .Machine$double.eps,
                  use_cpp = FALSE,
                  n_threads = 1) {
  lreplace(
    mmds(use_cpp = use_cpp, n_threads = n_threads),
    init = function(cost,
                    X,
                    max_iter,
                    eps = .Machine$double.eps,
                    verbose = FALSE,
                    ret_extra = c()) {
      cost <- mmds_init(
        cost,
        X,
        max_iter,
        eps,
        verbose,
        ret_extra,
        use_cpp = use_cpp,
        n_threads = n_threads
      )
      cost$R2 <- cost$R * cost$R
      cost$R <- NULL
      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)
      cost$pcost <- colSums((cost$R2 - cost$D2)^2)
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      cost$G <- k2g(Y, -8 * (cost$R2 - cost$D2))
      cost
    },
    export = function(cost, val) {
      res <- NULL
      switch(val,
        dx = {
          res <- cost$R2
          res[res < 0] <- 0
          sqrt(res)
        },
        dy = {
          res <- cost$D2
          res[res < 0] <- 0
          sqrt(res)
        }
      )
      res
    },
    update = function(cost, Y) {
      cost$D2 <- calc_d2(Y, use_cpp = use_cpp, n_threads = n_threads)
      cost
    },
    sentinel = "D2"
  )
}

sammon <- function(eps = .Machine$double.eps,
                   use_cpp = FALSE,
                   n_threads = 1) {
  lreplace(
    mmds(
      eps = eps,
      use_cpp = use_cpp,
      n_threads = n_threads
    ),
    init = function(cost,
                    X,
                    max_iter,
                    verbose = FALSE,
                    ret_extra = c()) {
      cost <- mmds_init(
        cost,
        X,
        max_iter,
        eps,
        verbose,
        ret_extra,
        use_cpp = use_cpp,
        n_threads = n_threads
      )
      cost$rsum_inv <- 1 / sum(cost$R)
      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)
      cost$pcost <- colSums((cost$R - cost$D)^2 / (cost$R + cost$eps)) * cost$rsum_inv
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      cost$G <- k2g(Y, -4 * cost$rsum_inv * (cost$R - cost$D) / (cost$R * cost$D + cost$eps))
      cost
    }
  )
}


gmmds <- function(k,
                  eps = .Machine$double.eps,
                  use_cpp = FALSE,
                  n_threads = 0) {
  lreplace(
    mmds(use_cpp = use_cpp, n_threads = n_threads),
    init = function(cost,
                    X,
                    max_iter,
                    verbose = FALSE,
                    ret_extra = c()) {
      cost$R <- geodesic(X,
        k,
        n_threads = n_threads,
        use_cpp = use_cpp,
        verbose = verbose
      )
      cost$eps <- eps
      cost
    },
    export = function(cost, val) {
      res <- NULL
      switch(val,
        geo = {
          res <- cost$R
        },
        dy = {
          res <- cost$D
        }
      )
      res
    }
  )
}

# Define neighborhoods using a radius based on a fraction (f) of all input
# distances (sorted by increasing length), don't correct non-neighborhood
# distances unless they smaller than the input distance
ballmmds <- function(f = 0.1,
                     eps = .Machine$double.eps,
                     use_cpp = FALSE,
                     n_threads = 1) {
  lreplace(
    mmds(use_cpp = use_cpp, n_threads = n_threads),
    init = function(cost,
                    X,
                    max_iter,
                    verbose = FALSE,
                    ret_extra = c()) {
      cost <- mmds_init(
        cost = cost,
        X = X,
        max_iter = max_iter,
        eps = eps,
        verbose = verbose,
        ret_extra = ret_extra,
        use_cpp = use_cpp,
        n_threads = n_threads
      )
      rs <- cost$R[upper.tri(cost$R)]
      rmax <- Rfast::nth(rs, max(1, round(f * length(rs))))
      if (verbose) {
        tsmessage("f = ", formatC(f), " rmax = ", formatC(rmax))
      }
      cost$rmax <- rmax

      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)

      R <- cost$R
      D <- cost$D
      rmax <- cost$rmax
      Ddiff <- R - D
      Ddiff[R > rmax & D > R] <- 0
      cost$pcost <- colSums(Ddiff * Ddiff)
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      eps <- cost$eps
      R <- cost$R
      D <- cost$D

      K <- -4 * (R - D) / (D + eps)

      rmax <- cost$rmax
      K[R > rmax & D > R] <- 0

      cost$G <- k2g(Y, K)
      cost
    },
    export = function(cost, val) {
      res <- NULL
      switch(val,
        geo = {
          res <- cost$R
        },
        dy = {
          res <- cost$D
        }
      )
      res
    },
    update = function(cost, Y) {
      cost$D <- calc_d(Y)
      cost
    }
  )
}

# Create the symmetrized knn graph, don't correct non-neighborhood distances
# unless they are smaller than the input distance
knnmmds <- function(k,
                    eps = .Machine$double.eps,
                    use_cpp = FALSE,
                    n_threads = 0) {
  lreplace(
    mmds(use_cpp = use_cpp, n_threads = n_threads),
    init = function(cost,
                    X,
                    max_iter,
                    verbose = FALSE,
                    ret_extra = c()) {
      cost <- mmds_init(
        cost = cost,
        X = X,
        max_iter = max_iter,
        eps = eps,
        verbose = verbose,
        ret_extra = ret_extra,
        use_cpp = use_cpp,
        n_threads = n_threads
      )
      knn <- knn_graph(
        X = X,
        k = k,
        ret_sparse = FALSE,
        n_threads = n_threads,
        verbose = verbose
      )
      # symmetrize
      cost$knn <- pmax(knn, t(knn))
      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)

      R <- cost$R
      D <- cost$D
      knn <- cost$knn
      Ddiff <- R - D
      Ddiff[knn == 0 & D > R] <- 0
      cost$pcost <- colSums(Ddiff * Ddiff)
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)

      eps <- cost$eps
      R <- cost$R
      D <- cost$D

      K <- -4 * (R - D) / (D + eps)

      knn <- cost$knn
      K[knn == 0 & D > R] <- 0

      cost$G <- k2g(Y, K)
      cost$D <- D
      cost
    },
    export = function(cost, val) {
      res <- NULL
      switch(val,
        geo = {
          res <- cost$R
        },
        dy = {
          res <- cost$D
        }
      )
      res
    },
    update = function(cost, Y) {
      cost$D <- calc_d(Y, use_cpp = use_cpp, n_threads = n_threads)
      cost
    }
  )
}

# Given data X and k nearest neighbors, return a geodisic distance matrix
# Disconnections are treated by using the Euclidean distance.
geodesic <- function(X,
                     k,
                     fill = TRUE,
                     use_cpp = FALSE,
                     n_threads = 0,
                     verbose = FALSE) {
  tsmessage("Calculating geodesic distances with k = ", k)

  R <- knn_dist(X, k, n_threads = n_threads, verbose = verbose)
  # The hard work is done by Rfast's implementation of Floyd's algorithm
  G <- Rfast::floyd(R)
  if (any(is.infinite(G)) && fill) {
    tsmessage(
      "k = ",
      k,
      " resulted in disconnections: filling with Euclidean distances"
    )
    if (methods::is(X, "dist")) {
      R <- as.matrix(X)
    } else {
      R <- calc_d(X, use_cpp = use_cpp, n_threads = n_threads)
    }
    G[is.infinite(G)] <- R[is.infinite(G)]
  }
  G
}
