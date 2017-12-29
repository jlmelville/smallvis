# Generic Functions -------------------------------------------------------

# Convert Force constant to Gradient
k2g <- function(Y, K, symmetrize = FALSE) {
  if (symmetrize) {
    K <- K + t(K)
  }
  Y * rowSums(K) - (K %*% Y)
}

cost_init <- function(cost, X, verbose = FALSE, ret_extra = c()) {
  if (!is.null(cost$init)) {
    cost <- cost$init(cost, X, verbose = verbose, ret_extra = ret_extra)
  }
  cost
}

cost_grad <- function(cost, Y) {
  cost$gr(cost, Y)
}

cost_point <- function(cost, Y) {
  cost$pfn(cost, Y)
}


# Cost Functions ----------------------------------------------------------

# LargeVis
# NB This version doesn't normalize the input P, despite what the paper
# indicates (source code of the current implementation doesn't seem to either)
largevis <- function(perplexity, gamma = 7, gr_eps = 0.1) {
  lreplace(tsne(perplexity),
     init = function(cost, X, eps = 1e-9, verbose = FALSE, ret_extra = c()) {
       cost <- sne_init(cost, X = X, perplexity = perplexity, symmetrize = "symmetric",
                        normalize = FALSE, verbose = verbose,
                        ret_extra = ret_extra)
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
     }
  )
}

# UMAP
umap <- function(perplexity, spread = 1, min_dist = 0.001, gr_eps = 0.001) {
  lreplace(tsne(perplexity),
    init = function(cost, X, eps = 1e-9, verbose = FALSE, ret_extra = c()) {

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
    }
  )
}

# UMAP with the output kernel fixed to the t-distribution
tumap <- function(perplexity, gr_eps = 0.1) {
  lreplace(tsne(perplexity),
    init = function(cost, X, eps = 1e-9, verbose = FALSE, ret_extra = c()) {
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
    }
  )
}

# t-UMAP where output and input affinities are normalized
ntumap <- function(perplexity, gr_eps = 0.1) {
  lreplace(tsne(perplexity),
    init = function(cost, X, eps = 1e-9, verbose = FALSE, ret_extra = c()) {

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
    }
  )
}

mmds_init <- function(cost, X, eps = .Machine$double.eps, verbose = FALSE,
                      ret_extra = c()) {
  if (methods::is(X, "dist")) {
    cost$R <- X
  }
  else {
    cost$R <- sqrt(safe_dist2(X))
  }
  cost$eps <- eps
  cost
}

# Metric MDS, minimizing strain.
mmds <- function() {
  list(
    init = mmds_init,
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

smmds <- function() {
  lreplace(
    mmds(),
    init = function(cost, X, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      cost <- mmds_init(cost, X, eps, verbose, ret_extra)
      cost$R2 <- cost$R * cost$R
      cost$R <- NULL
      cost
    },
    pfn = function(cost, Y) {
      cost$pcost <- colSums((cost$R2 - cost$D2) ^ 2)
      cost
    },
    gr = function(cost, Y) {
      D2 <- dist2(Y)
      cost$G <- k2g(Y,  -8 * (cost$R2 - D2))
      cost$D2 <- D2
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
    })
}

sammon <- function() {
  lreplace(mmds(),
    init = function(cost, X, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      cost <- mmds_init(cost, X, eps, verbose, ret_extra)
      cost$rsum_inv <- 1 / sum(cost$R)
      cost
    },
    pfn = function(cost, Y) {
      eps <- cost$eps
      R <- cost$R
      D <- cost$D
      cost$pcost <- colSums((R - D) ^ 2 / (R + eps)) * cost$rsum_inv
      cost
    },
    gr = function(cost, Y) {
      eps <- cost$eps
      R <- cost$R
      D <- sqrt(safe_dist2(Y))
      cost$G <- k2g(Y,  -4 * cost$rsum_inv * (R - D) / (R * D + eps))
      cost$D <- D
      cost
    }
  )
}


geommds <- function(k) {
  lreplace(
    mmds(),
    init = function(cost, X, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
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

# Carreira-Perpinán, M. A. (2010, June).
# The Elastic Embedding Algorithm for Dimensionality Reduction.
# In \emph{Proceedings of the 27th International Conference on Machine Learning (ICML-10)} (pp. 167-174).
# http://faculty.ucmerced.edu/mcarreira-perpinan/papers/icml10.pdf (PDF)
ee <- function(perplexity, lambda = 100) {
  list(
    init = function(cost, X, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      if (methods::is(X, "dist")) {
        R <- X
      }
      else {
        R <- sqrt(safe_dist2(X))
      }
      cost$Vn <- R / sum(R)
      cost <- sne_init(cost, X, perplexity = perplexity,
                       symmetrize = "symmetric", normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra)
      cost$eps <- eps
      cost
    },
    pfn = function(cost, Y) {
      Vp <- cost$P
      Vn <- cost$Vn
      W <- cost$W
      eps <- cost$eps
      cost$pcost <- colSums(-Vp * log(W + eps) + lambda * (Vn * W))
      cost
    },
    gr = function(cost, Y) {
      Vp <- cost$P
      Vn <- cost$Vn

      W <- dist2(Y)
      W <- exp(-W)
      diag(W) <- 0
      cost$W <- W
      cost$G <- k2g(Y,  4 * (Vp - lambda * Vn * W))
      cost
    },
    export = function(cost, val) {
      res <- NULL
      if (!is.null(cost[[val]])) {
        res <- cost[[val]]
      }
      else if (!is.null(cost[[toupper(val)]])) {
        res <- cost[[toupper(val)]]
      }
      res
    }
  )
}


# JSE and NeRV ------------------------------------------------------------

# Venna, J., Peltonen, J., Nybo, K., Aidos, H., & Kaski, S. (2010).
# Information retrieval perspective to nonlinear dimensionality reduction for
# data visualization.
# \emph{Journal of Machine Learning Research}, \emph{11}, 451-490.
#
# Unlike original publication, won't transfer input precisions to output kernel
# lambda = 1 gives ASNE results
# default lambda = 0.9 from "Majorization-Minimization for Manifold Embedding"
# Yang, Peltonen, Kaski 2015
nerv <- function(perplexity, lambda = 0.9) {
  lreplace(
    tsne(perplexity = perplexity),
    init = function(cost, X, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      cost <- sne_init(cost, X, perplexity = perplexity,
                       symmetrize = "none", normalize = FALSE,
                       verbose = verbose, ret_extra = ret_extra)
      cost$eps <- eps
      cost
    },
    pfn = function(cost, Y) {
      P <- cost$P

      kl_fwd <- rowSums(P * cost$lPQ)

      cost$pcost <- lambda * kl_fwd + (1 - lambda) * cost$kl_rev
      cost
    },
    gr = function(cost, Y) {
      eps <- cost$eps
      P <- cost$P
      W <- dist2(Y)
      W <- exp(-W)

      # Particularly for low lambda, need to make sure W is never all zero
      W[W < eps] <- eps
      diag(W) <- 0
      invZ <- 1 / colSums(W)
      Q <- W * invZ

      # Forward KL gradient
      K <- lambda * (P - Q)

      # Reverse KL gradient
      lPQ <- log((P + eps) / (Q + eps))
      # for KLrev we want Q * log(Q/P), so take -ve of log(P/Q)
      kl_rev <- rowSums(Q * -lPQ)

      # Total K
      K <- K + (1 - lambda) * (Q * (lPQ + kl_rev))

      cost$G <- k2g(Y, 2 * K, symmetrize = TRUE)
      cost$invZ <- invZ
      cost$W <- W
      cost$kl_rev <- kl_rev
      cost$lPQ <- lPQ

      cost
    })
}

# Lee, J. A., Renard, E., Bernard, G., Dupont, P., & Verleysen, M. (2013).
# Type 1 and 2 mixtures of Kullback-Leibler divergences as cost functions in
# dimensionality reduction based on similarity preservation.
# \emph{Neurocomputing}, \emph{112}, 92-108.
# kappa = 0 behaves like ASNE
# kappa = 1 behaves like NeRV with lambda = 0. Yes that's confusing.
jse <- function(perplexity, kappa = 0.5) {
  # safeguard against 0 kappa, because gradient requires kappa to be > 0
  # check for kappa == 1 for better accuracy
  kappa <- max(1e-8, kappa)
  kappa_inv <- 1 / kappa
  om_kappa <- 1 - kappa
  om_kappa_inv <- 1 / om_kappa

  lreplace(
    tsne(perplexity = perplexity),
    init = function(cost, X, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      cost <- sne_init(cost, X, perplexity = perplexity,
                       symmetrize = "none", normalize = FALSE,
                       verbose = verbose, ret_extra = ret_extra)
      cost$eps <- eps
      cost
    },
    pfn = function(cost, Y) {
      eps <- cost$eps
      P <- cost$P
      Z <- cost$Z

      kl_fwd <- rowSums(P * log((P + eps) /  (Z + eps)))

      if (kappa == 1) {
        cost$pcost <- cost$kl_rev
      }
      else {
        cost$pcost <- om_kappa_inv * kl_fwd + kappa_inv * cost$kl_rev
      }

      cost
    },
    gr = function(cost, Y) {
      eps <- cost$eps
      P <- cost$P
      W <- dist2(Y)
      W <- exp(-W)

      W[W < eps] <- eps
      diag(W) <- 0
      # nomenclature overlap problem here
      # normally sum of W is Z, but in JSE that's the combination of P and Q
      invS <- 1 / colSums(W)
      Q <- W * invS

      if (kappa == 1) {
        Z <- P
      }
      else {
        Z <- kappa * P + om_kappa * Q
      }

      lZQ <- log((Z + eps) / (Q + eps))
      kl_rev <- rowSums(Q * -lZQ)
      K <- kappa_inv * (Q * (lZQ + kl_rev))

      cost$G <- k2g(Y, 2 * K, symmetrize = TRUE)
      cost$invZ <- invS
      cost$W <- W
      cost$kl_rev <- kl_rev
      cost$Z <- Z

      cost
    })
}


