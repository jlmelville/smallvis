# Cost Functions ----------------------------------------------------------

# UMAP
umap <- function(perplexity, spread = 1, min_dist = 0.001, gr_eps = 0.001) {
  list(
    init = function(cost, X, eps = 1e-9, verbose = FALSE) {

      ab_params <- find_ab_params(spread = spread, min_dist = min_dist)
      a <- ab_params[1]
      b <- ab_params[2]
      if (verbose) {
        message("Umap curve parameters = ", formatC(a), ", ", formatC(b))
      }

      P <- smooth_knn_distances(X, k = perplexity, tol = 1e-5,
                                verbose = verbose)$P
      # Fuzzy set union
      P <- P + t(P) - P * t(P)
      cost$P <- P

      cost$a <- a
      cost$b <- b
      cost$eps <- eps

      cost
    },
    pfn = function(cost, Y) {
      P <- cost$P
      eps <- cost$eps
      W <- cost$W
      cost$pcost <- colSums(-P * log(W + eps) - (1 - P) * log1p(-W + eps))
      cost
    },
    gr = function(cost, Y) {
      P <- cost$P
      a <- cost$a
      b <- cost$b
      eps <- cost$eps

      D2 <- dist2(Y)
      D2[D2 < 0] <- 0

      W <- 1 / (1 + a * D2 ^ b)
      diag(W) <- 0

      WF <- a * b * (D2 + eps) ^ (b - 1)
      diag(WF) <- 0
      WF <- W * WF

      cost$G <- k2g(Y, 4 * (P * WF - ((1 - P) * W * WF) / ((1 - W) + gr_eps)))
      cost$W <- W
      cost
    },
    export = function(cost, val) {
      res <- NULL
      switch(val,
             w = {
               res <- cost$W
             },
             p = {
               res <- cost$P
             }
      )
      res
    }
  )
}

# LargeVis
# NB This version doesn't normalize the input P, despite what the paper
# indicates (source code of the current implementation doesn't seem to either)
largevis <- function(perplexity, gamma = 7, gr_eps = 0.1) {
  list(
    init = function(cost, X, eps = 1e-9, verbose = FALSE) {
      P <- sne_init(X = X, perplexity = perplexity, symmetrize = "symmetric",
                    normalize = FALSE, verbose = verbose)
      cost$eps <- eps
      cost
    },
    pfn = function(cost, Y) {
      P <- cost$P
      eps <- cost$eps
      W <- cost$W
      cost$pcost <- colSums(-P * log(W + eps) - gamma * log1p(-W + eps))
      cost
    },
    gr = function(cost, Y) {
      P <- cost$P
      W <- dist2(Y)

      W <- 1 / (1 + W)
      diag(W) <- 0
      cost$G <- k2g(Y, 4 * (P * W - ((gamma * W * W) / ((1 - W) + gr_eps))))
      cost$W <- W
      cost
    },
    export = function(cost, val) {
      res <- NULL
      switch(val,
             w = {
               res <- cost$W
             },
             p = {
               res <- cost$P
             }
      )
      res
    }
  )
}

# UMAP with the output kernel fixed to the t-distribution
tumap <- function(perplexity, gr_eps = 0.1) {
  list(
    init = function(cost, X, eps = 1e-9, verbose = FALSE) {
      cost$eps <- eps

      P <- smooth_knn_distances(X, k = perplexity, tol = 1e-5,
                                verbose = verbose)$P
      # Fuzzy set union
      P <- P + t(P) - P * t(P)
      cost$P <- P

      cost
    },
    pfn = function(cost, Y) {
      P <- cost$P
      eps <- cost$eps
      W <- cost$W
      cost$pcost <- colSums(-P * log(W + eps) - (1 - P) * log1p(-W + eps))
      cost
    },
    gr = function(cost, Y) {
      P <- cost$P
      W <- dist2(Y)
      W <- 1 / (1 + W)
      diag(W) <- 0

      cost$G <- k2g(Y,  4 * (P * W - ((1 - P) * W * W) / ((1 - W) + gr_eps)))
      cost$W <- W
      cost
    },
    export = function(cost, val) {
      res <- NULL
      switch(val,
             w = {
               res <- cost$W
             },
             p = {
               res <- cost$P
             }
      )
      res
    }
  )
}

# t-UMAP where output and input affinities are normalized
ntumap <- function(perplexity, gr_eps = 0.1) {
  list(
    init = function(cost, X, eps = 1e-9, verbose = FALSE) {

      P <- smooth_knn_distances(X, k = perplexity, tol = 1e-5,
                                verbose = verbose)$P
      # Symmetrize by fuzzy set union
      P <- P + t(P) - P * t(P)
      # Normalize
      cost$P <- P / sum(P)
      cost$eps <- eps
      cost
    },
    pfn = function(cost, Y) {
      P <- cost$P
      eps <- cost$eps
      W <- cost$W
      Q <- W * cost$invZ

      cost$pcost <- colSums(-P * log(Q + eps) - (1 - P) * log1p(-Q + eps))
      cost
    },
    gr = function(cost, Y) {
      P <- cost$P
      W <- dist2(Y)
      W <- 1 / (1 + W)
      diag(W) <- 0

      invZ <- 1 / sum(W)
      Q <- W * invZ
      C <- (P - Q) / (1 - Q)

      cost$G <- k2g(Y,  4 * W * (C - sum(C) * Q))
      cost$W <- W
      cost$invZ <- invZ
      cost
    },
    export = function(cost, val) {
      res <- NULL
      switch(val,
             w = {
               res <- cost$W
             },
             q = {
               res <- cost$W * cost$invZ
             },
             p = {
               res <- cost$P
             }
      )
      res
    }
  )
}

# Metric MDS, minimizing strain.
mmds <- function() {
  list(
    init = function(cost, X, eps = .Machine$double.eps, verbose = FALSE) {
      if (methods::is(X, "dist")) {
        cost$R <- X
      }
      else {
        cost$R <- sqrt(safe_dist2(X))
      }
      cost$eps <- eps
      cost
    },
    pfn = function(cost, Y) {
      R <- cost$R
      D <- cost$D
      cost$pcost <- colSums((R - D) ^ 2)
      cost
    },
    gr = function(cost, Y) {
      eps <- cost$eps
      R <- cost$R
      D <- sqrt(safe_dist2(Y))
      cost$G <- k2g(Y,  -4 * (R - D) / (D + eps))
      cost$D <- D
      cost
    },
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

geommds <- function(k) {
  lreplace(
    mmds(),
    init = function(cost, X, eps = .Machine$double.eps, verbose = FALSE) {
      cost$R <- geodesic(X, k)

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

# Generic Functions -------------------------------------------------------

# Convert Force constant to Gradient
k2g <- function(Y, K, symmetrize = FALSE) {
  if (symmetrize) {
    K <- K + t(K)
  }
  Y * rowSums(K) - (K %*% Y)
}

cost_init <- function(cost, X, verbose = FALSE) {
  if (!is.null(cost$init)) {
    cost <- cost$init(cost, X, verbose = verbose)
  }
  cost
}

cost_grad <- function(cost, Y) {
  cost$gr(cost, Y)
}

cost_point <- function(cost, Y) {
  cost$pfn(cost, Y)
}
