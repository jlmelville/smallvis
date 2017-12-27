# t-SNE
tsne <- function(perplexity, inp_kernel = "gaussian") {
  list(
    init = function(cost, X, eps = .Machine$double.eps, verbose = FALSE) {
      cost$P <- sne_init(X, perplexity = perplexity, inp_kernel = inp_kernel,
                         symmetrize = "symmetric", normalize = TRUE,
                         verbose = verbose)
      cost$eps <- eps
      cost
    },
    pfn = function(cost, Y) {
      P <- cost$P
      eps <- cost$eps
      invZ <- cost$invZ
      W <- cost$W
      cost$pcost <- colSums(P * log((P + eps) / ((W * invZ) + eps)))
      cost
    },
    gr = function(cost, Y) {
      P <- cost$P
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
             },
             p = {
               res <- cost$P
             }
      )
      res
    }
  )
}

# Cook, J., Sutskever, I., Mnih, A., & Hinton, G. E. (2007).
# Visualizing similarity data with a mixture of maps.
# In \emph{International Conference on Artificial Intelligence and Statistics} (pp. 67-74).
ssne <- function(perplexity, inp_kernel = "gaussian") {
  lreplace(
    tsne(perplexity = perplexity, inp_kernel = inp_kernel),
    gr = function(cost, Y) {
      P <- cost$P
      W <- dist2(Y)
      W <- exp(-W)
      diag(W) <- 0
      invZ <- 1 / sum(W)
      cost$invZ <- invZ
      cost$W <- W
      cost$G <- k2g(Y, 4 * (P - W * invZ))
      cost
    }
  )
}

# Hinton, G. E., & Roweis, S. T. (2002).
# Stochastic neighbor embedding.
# In \emph{Advances in neural information processing systems} (pp. 833-840).
asne <- function(perplexity, inp_kernel = "gaussian") {
  list(
    init = function(cost, X, eps = .Machine$double.eps, verbose = FALSE) {
      cost$P <- sne_init(X, perplexity = perplexity, inp_kernel = inp_kernel,
                         symmetrize = "none", normalize = FALSE,
                         verbose = verbose)
      cost$eps <- eps
      cost
    },
    pfn = function(cost, Y) {
      P <- cost$P
      eps <- cost$eps
      invZ <- cost$invZ
      W <- cost$W

      # ASNE defines N KL divergences, row-wise
      cost$pcost <- rowSums(P * log((P + eps) / ((W * invZ) + eps)))
      cost
    },
    gr = function(cost, Y) {
      P <- cost$P
      W <- dist2(Y)
      W <- exp(-W)
      diag(W) <- 0
      invZ <- 1 / colSums(W)
      K <- (P - W * invZ)
      cost$G <- k2g(Y, 2 * K, symmetrize = TRUE)
      cost$invZ <- invZ
      cost$W <- W

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

# Heavy-Tailed Symmetric Stochastic Neighbor Embedding (HSSNE)
# Yang, Z., King, I., Xu, Z., & Oja, E. (2009).
# Heavy-tailed symmetric stochastic neighbor embedding.
# In \emph{Advances in neural information processing systems} (pp. 2169-2177).
hssne <- function(perplexity, inp_kernel = "gaussian", alpha = 0.5) {
  lreplace(
    tsne(perplexity = perplexity, inp_kernel = inp_kernel),
    gr = function(cost, Y) {
      P <- cost$P
      W <- dist2(Y)
      # to include bandwidth
      # W <- (alpha * beta * W + 1) ^ (-1 / alpha)
      W <- (alpha * W + 1) ^ (-1 / alpha)
      diag(W) <- 0

      invZ <- 1 / sum(W)
      cost$invZ <- invZ
      cost$W <- W
      # to include bandwidth
      # K <- 4 * beta * (P - W * invZ) * (W ^ alpha)
      cost$G <- k2g(Y, 4 * (P - W * invZ) * (W ^ alpha))
      cost
    }
  )
}

wtsne <- function(perplexity, inp_kernel = "gaussian") {
  list(
    init = function(cost, X, eps = .Machine$double.eps, verbose = FALSE) {
      cost$P <- sne_init(X, perplexity = perplexity, inp_kernel = inp_kernel,
                         symmetrize = "symmetric", normalize = TRUE,
                         verbose = verbose)
      # degree centrality
      deg <- colSums(cost$P)
      cost$M <- outer(deg, deg) * nrow(cost$P)
      cost$eps <- eps
      cost
    },
    pfn = function(cost, Y) {
      P <- cost$P
      eps <- cost$eps
      invZ <- cost$invZ
      W <- cost$W
      cost$pcost <- colSums(P * log((P + eps) / ((W * invZ) + eps)))
      cost
    },
    gr = function(cost, Y) {
      P <- cost$P
      M <- cost$M
      W <- dist2(Y)
      W <- M / (1 + W)
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
             },
             p = {
               res <- cost$P
             }
      )
      res
    }
  )
}

wssne <- function(perplexity) {
  lreplace(
    wtsne(perplexity = perplexity),
    gr = function(cost, Y) {
      P <- cost$P
      M <- cost$M
      W <- dist2(Y)
      W <- M * exp(-W)
      diag(W) <- 0

      invZ <- 1 / sum(W)
      cost$invZ <- invZ
      cost$W <- W
      cost$G <- k2g(Y, 4 * (P - W * invZ))
      cost
    }
  )
}


# symmetrize: type of symmetrization:
#  none - no symmetrization as in ASNE, JSE, NeRV
#  symmetric - symmetric nearest neighbor style, default, as in t-SNE.
#  mutual - mutual nearest neighbor style as suggested by Schubert and Gertz in
#  "Intrinsic t-Stochastic Neighbor Embedding for Visualization and Outlier
#   Detection - A Remedy Against the Curse of Dimensionality?"
sne_init <- function(X, perplexity, inp_kernel = "gaussian",
                     symmetrize = "symmetric",
                     normalize = TRUE,
                     verbose = FALSE) {
  P <- x2p(X, perplexity, tol = 1e-5, kernel = inp_kernel, verbose = verbose)$P

  # Symmetrize
  P <- switch(symmetrize,
              none = P,
              symmetric = 0.5 * (P + t(P)),
              mutual = sqrt(P * t(P)),
              stop("unknown symmetrization: ", symmetrize))
  # Normalize
  if (normalize) {
    P <- P / sum(P)
  }
  P
}
