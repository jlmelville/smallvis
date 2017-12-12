# Cost Functions ----------------------------------------------------------

# t-SNE
tsne <- function() {
  list(
    init = function(cost, X, eps = .Machine$double.eps, verbose = FALSE) {
      cost$eps <- eps
      cost
    },
    pfn = function(cost, P, Y) {
      eps <- cost$eps
      invZ <- cost$invZ
      W <- cost$W
      cost$pcost <- colSums(P * log((P + eps) / ((W * invZ) + eps)))
      cost
    },
    gr = function(cost, P, Y) {
      W <- dist2(Y)
      W <- 1 / (1 + W)
      diag(W) <- 0
      invZ <- 1 / sum(W)
      cost$invZ <- invZ
      cost$W <- W
      cost$G <- k2g(Y, 4 * W * (P - W * invZ))
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
        }
      )
      res
    }
  )
}

# UMAP
umap <- function(spread = 1, min_dist = 0.001, gr_eps = 0.1) {
  list(
    init = function(cost, X, eps = 1e-9, verbose = FALSE) {

      ab_params <- find_ab_params(spread = spread, min_dist = min_dist)
      a <- ab_params[1]
      b <- ab_params[2]
      if (verbose) {
        message("Umap curve parameters = ", formatC(a), ", ", formatC(b))
      }

      cost$a <- a
      cost$b <- b
      cost$eps <- eps

      cost
    },
    pfn = function(cost, P, Y) {
      eps <- cost$eps
      W <- cost$W
      cost$pcost <- colSums(-P * log(W + eps) - (1 - P) * log1p(-W + eps))
      cost
    },
    gr = function(cost, P, Y) {
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
             }
      )
      res
    }
  )
}

# LargeVis
largevis <- function(gamma = 7, gr_eps = 0.1) {
  list(
    init = function(cost, X, eps = 1e-9, verbose = FALSE) {
      cost$eps <- eps
      cost
    },
    pfn = function(cost, P, Y) {
      eps <- cost$eps
      W <- cost$W
      cost$pcost <- colSums(-P * log(W + eps) - gamma * log1p(-W + eps))
      cost
    },
    gr = function(cost, P, Y) {
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
             }
      )
      res
    }
  )
}

# UMAP with the output kernel fixed to the t-distribution
tumap <- function(gr_eps = 0.1) {
  list(
    init = function(cost, X, eps = 1e-9, verbose = FALSE) {
      cost$eps <- eps
      cost
    },
    pfn = function(cost, P, Y) {
      eps <- cost$eps
      W <- cost$W
      cost$pcost <- colSums(-P * log(W + eps) - (1 - P) * log1p(-W + eps))
      cost
    },
    gr = function(cost, P, Y) {
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
             }
      )
      res
    }
  )
}

# t-UMAP where output and input affinities are normalized
ntumap <- function(gr_eps = 0.1) {
  list(
    init = function(cost, X, eps = 1e-9, verbose = FALSE) {
      cost$eps <- eps
      cost
    },
    pfn = function(cost, P, Y) {
      eps <- cost$eps
      W <- cost$W
      Q <- W * cost$invZ

      cost$pcost <- colSums(-P * log(Q + eps) - (1 - P) * log1p(-Q + eps))
      cost
    },
    gr = function(cost, P, Y) {
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
             }
      )
      res
    }
  )
}

# Generic Functions -------------------------------------------------------

# Convert Force constant to Gradient
k2g <- function(Y, K) {
  Y * rowSums(K) - (K %*% Y)
}

cost_init <- function(cost, X, verbose = FALSE) {
  if (!is.null(cost$init)) {
    cost <- cost$init(cost, X, verbose = verbose)
  }
  cost
}

cost_grad <- function(cost, P, Y) {
  cost$gr(cost, P, Y)
}

cost_point <- function(cost, P, Y) {
  cost$pfn(cost, P, Y)
}
