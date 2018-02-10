kl_cost <- function(cost, Y) {
  P <- cost$P
  eps <- cost$eps
  invZ <- cost$invZ
  W <- cost$W
  # P log(P / Q) = P log P - P log Q
  cost$pcost <- cost$plogp - colSums(P * log(((W * invZ) + eps)))
  cost
}

# t-SNE
tsne <- function(perplexity, inp_kernel = "gaussian") {
  list(
    init = function(cost, X, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = "symmetric", normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra)
      P <- cost$P
      # cache P log P constant part of cost: incur 1 extra log operation now
      # but saves one division operation every time we calculate cost:
      # substantial (10-15%) speed up with Wolfe line search methods
      cost$plogp <- colSums(P * log((P + eps)))

      cost$eps <- eps
      cost
    },
    pfn = kl_cost,
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
      if (!is.null(cost[[val]])) {
        res <- cost[[val]]
      }
      else if (!is.null(cost[[toupper(val)]])) {
        res <- cost[[toupper(val)]]
      }
      else {
        switch(val,
               q = {
                 res <- cost$W * cost$invZ
               })
      }
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
  lreplace(tsne(perplexity),
    init = function(cost, X, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = "none", normalize = FALSE,
                       verbose = verbose, ret_extra = ret_extra)
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
      cost$G <- k2g(Y, 2 * (P - W * invZ), symmetrize = TRUE)
      cost$invZ <- invZ
      cost$W <- W

      cost
    }
  )
}

# Heavy-Tailed Symmetric Stochastic Neighbor Embedding (HSSNE)
# Yang, Z., King, I., Xu, Z., & Oja, E. (2009).
# Heavy-tailed symmetric stochastic neighbor embedding.
# In \emph{Advances in neural information processing systems} (pp. 2169-2177).
hssne <- function(perplexity, alpha = 0.5) {
  alpha <- max(alpha, 1e-8)
  lreplace(
    tsne(perplexity = perplexity),
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

# A version of HSSNE where alpha is allowed to vary at every epoch
dhssne <- function(perplexity, alpha = 0.5) {
  alpha_min <- 1e-8
  alpha <- max(alpha, alpha_min)
  lreplace(
    tsne(perplexity = perplexity),
    gr = function(cost, Y) {
      alpha <- cost$alpha
      P <- cost$P
      W <- dist2(Y)
      W <- (alpha * W + 1) ^ (-1 / alpha)
      diag(W) <- 0

      invZ <- 1 / sum(W)
      cost$invZ <- invZ
      cost$W <- W

      cost$G <- k2g(Y, 4 * (P - W * invZ) * (W ^ alpha))
      cost
    },
    epoch = function(opt, cost, iter, Y) {
      errs <- c()

      old_cost <- cost

      alpha_curr <- cost$alpha

      alpha_down <- max(alpha_curr - 0.1, alpha_min)
      alpha_up <- alpha_curr + 0.1
      alphas <- unique(c(alpha_down, alpha_curr, alpha_up))

      for (alpha in alphas) {
        cost$alpha <- alpha
        cost <- cost$gr(cost, Y)
        cost <- cost_point(cost, Y)
        err <- sum(cost$pcost)
        errs <- c(errs, err)
      }

      old_cost$alpha <- alphas[which.min(errs)]
      list(cost = old_cost)
    },
    alpha = alpha
  )
}

# Yang, Z., Peltonen, J., & Kaski, S. (2014).
# Optimization equivalence of divergences improves neighbor embedding.
# In \emph{Proceedings of the 31st International Conference on Machine Learning (ICML-14)}
# (pp. 460-468).
wtsne <- function(perplexity) {
  lreplace(tsne(perplexity = perplexity),
    init = function(cost, X, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      ret_extra <- c(ret_extra, "pdeg")
      cost <- sne_init(cost, X, perplexity = perplexity,
                         symmetrize = "symmetric", normalize = TRUE,
                         verbose = verbose, ret_extra = ret_extra)
      # P matrix degree centrality: column sums
      deg <- cost$pdeg
      if (verbose) {
        summarize(deg, "deg")
      }
      cost$M <- outer(deg, deg)
      cost$invM <- 1 / cost$M
      cost$eps <- eps

      cost$plogp <- colSums(cost$P * log((cost$P + eps)))
      cost
    },
    gr = function(cost, Y) {
      P <- cost$P
      invM <- cost$invM
      M <- cost$M

      W <- dist2(Y)
      W <- M / (1 + W)
      diag(W) <- 0
      invZ <- 1 / sum(W)

      cost$invZ <- invZ
      cost$W <- W
      cost$G <- k2g(Y, 4 * W * invM * (P - W * invZ))
      cost
    }
  )
}

wssne <- function(perplexity) {
  lreplace(wtsne(perplexity = perplexity),
     gr = function(cost, Y) {
       P <- cost$P
       invM <- cost$invM
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


# Bandwidth Experiments --------------------------------------------------

# t-SNE with input kernel bandwidths transferred to output
btsne <- function(perplexity, inp_kernel = "gaussian") {
  lreplace(tsne(perplexity = perplexity, inp_kernel = inp_kernel),
    init = function(cost, X, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      ret_extra <- unique(c(ret_extra, 'beta'))
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = "symmetric", normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra)
      P <- cost$P
      cost$plogp <- colSums(P * log((P + eps)))

      cost$eps <- eps
      cost
    },
    gr = function(cost, Y) {
      P <- cost$P
      beta <- cost$beta
      W <- dist2(Y)
      W <- 1 / (1 + (beta * W))
      diag(W) <- 0
      invZ <- 1 / sum(W)
      cost$invZ <- invZ
      cost$W <- W
      cost$G <- k2g(Y, 2 * beta * W * (P - W * invZ), symmetrize = TRUE)
      cost
    }
  )
}

# SSNE with input kernel bandwidths transferred to output
bssne <- function(perplexity, inp_kernel = "gaussian") {
  lreplace(btsne(perplexity = perplexity, inp_kernel = inp_kernel),
           gr = function(cost, Y) {
             P <- cost$P
             beta <- cost$beta
             W <- dist2(Y)
             W <- exp(-beta * W)
             diag(W) <- 0
             invZ <- 1 / sum(W)
             cost$invZ <- invZ
             cost$W <- W
             cost$G <- k2g(Y, 2 * beta * (P - W * invZ), symmetrize = TRUE)
             cost
           }
  )
}

# ASNE with input kernel bandwidths transferred to output
basne <- function(perplexity) {
  lreplace(asne(perplexity = perplexity),
           init = function(cost, X, eps = .Machine$double.eps, verbose = FALSE,
                           ret_extra = c()) {
             ret_extra <- unique(c(ret_extra, 'beta'))

             cost <- sne_init(cost, X, perplexity = perplexity,
                              symmetrize = "none", normalize = FALSE,
                              verbose = verbose, ret_extra = ret_extra)
             cost$eps <- eps
             cost
           },
           gr = function(cost, Y) {
             P <- cost$P
             beta <- cost$beta
             W <- dist2(Y)
             W <- exp(-W * beta)
             diag(W) <- 0
             # NB Can't use colSums now that W is not symmetric!
             invZ <- 1 / rowSums(W)
             cost$G <- k2g(Y, 2 * beta * (P - W * invZ), symmetrize = TRUE)

             cost$invZ <- invZ
             cost$W <- W
             cost
           }
  )
}

# t-ASNE with input kernel bandwidths transferred to output
btasne <- function(perplexity) {
  lreplace(basne(perplexity = perplexity),
           gr = function(cost, Y) {
             beta <- cost$beta
             P <- cost$P
             W <- dist2(Y)
             W <- 1 / (1 + beta * W)
             diag(W) <- 0
             invZ <- 1 / rowSums(W)
             cost$invZ <- invZ
             cost$W <- W
             cost$G <- k2g(Y, 2 * beta * W * (P - W * invZ), symmetrize = TRUE)
             cost$invZ <- invZ
             cost$W <- W
             cost
           }
  )
}

# Normalization Experiments -----------------------------------------------

# ASNE but with the t-distributed kernel
tasne <- function(perplexity) {
  lreplace(asne(perplexity = perplexity),
  gr = function(cost, Y) {
    P <- cost$P
    W <- dist2(Y)
    W <- 1 / (1 + W)
    diag(W) <- 0
    invZ <- 1 / rowSums(W)
    cost$invZ <- invZ
    cost$W <- W
    cost$G <- k2g(Y, 2 * W * (P - W * invZ), symmetrize = TRUE)
    cost$invZ <- invZ
    cost$W <- W
    cost
  }
  )
}

# t-RM-SNE
# t-SNE without symmetrization of P (but still pair-normalizing)
# row-normalize, then matrix normalize
trmsne <- function(perplexity, inp_kernel = "gaussian") {
  lreplace(tsne(perplexity = perplexity, inp_kernel = inp_kernel),
           init = function(cost, X, eps = .Machine$double.eps, verbose = FALSE,
                           ret_extra = c()) {
             cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                              symmetrize = "none", normalize = TRUE,
                              verbose = verbose, ret_extra = ret_extra)
             P <- cost$P
             cost$plogp <- colSums(P * log((P + eps)))

             cost$eps <- eps
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
             cost$G <- k2g(Y, 2 * W * (P - W * invZ), symmetrize = TRUE)
             cost
           }
  )
}

# t-M-SNE
# t-SNE but without row-normalizing or symmetrizing, just matrix normalization
# Not recommended
tmsne <- function(perplexity, inp_kernel = "gaussian") {
  lreplace(trmsne(perplexity = perplexity, inp_kernel = inp_kernel),
           init = function(cost, X, eps = .Machine$double.eps, verbose = FALSE,
                           ret_extra = c()) {
             cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                              symmetrize = "none", row_normalize = FALSE,
                              normalize = TRUE,
                              verbose = verbose, ret_extra = ret_extra)
             P <- cost$P
             cost$plogp <- colSums(P * log((P + eps)))

             cost$eps <- eps
             cost
           }
  )
}

# RSR row-normalize, symmetrize, then row-normalize again
# Might work a tiny bit better than t-ASNE?
trsrsne <- function(perplexity) {
  lreplace(tasne(perplexity),
           init = function(cost, X, eps = .Machine$double.eps, verbose = FALSE,
                           ret_extra = c()) {
             cost <- sne_init(cost, X, perplexity = perplexity,
                              symmetrize = "symmetric", normalize = FALSE,
                              verbose = verbose, ret_extra = ret_extra)
             P <- cost$P
             P <- P / rowSums(P)
             cost$P <- P

             cost$eps <- eps
             cost
           }
  )
}


arsrsne <- function(perplexity) {
  lreplace(asne(perplexity),
           init = function(cost, X, eps = .Machine$double.eps, verbose = FALSE,
                           ret_extra = c()) {
             cost <- sne_init(cost, X, perplexity = perplexity,
                              symmetrize = "symmetric", normalize = FALSE,
                              verbose = verbose, ret_extra = ret_extra)
             P <- cost$P
             P <- P / rowSums(P)
             cost$P <- P

             cost$eps <- eps
             cost
           }
  )
}

# Perplexity Calibration --------------------------------------------------

# symmetrize: type of symmetrization:
#  none - no symmetrization as in ASNE, JSE, NeRV
#  symmetric - symmetric nearest neighbor style, default, as in t-SNE.
#  mutual - mutual nearest neighbor style as suggested by Schubert and Gertz in
#  "Intrinsic t-Stochastic Neighbor Embedding for Visualization and Outlier
#   Detection - A Remedy Against the Curse of Dimensionality?"
sne_init <- function(cost, X, perplexity, kernel = "gaussian",
                     symmetrize = "symmetric", row_normalize = TRUE,
                     normalize = TRUE,
                     verbose = FALSE, ret_extra = c()) {
  if (tolower(kernel) == "knn") {
    if (verbose) {
      tsmessage("Using knn kernel with k = ", formatC(perplexity))
    }
    P <- knn_graph(X, k = perplexity)
    x2ares <- list(W = P)
  }
  else {
    if (verbose) {
      tsmessage("Commencing calibration for perplexity = ", formatC(perplexity))
    }
    x2ares <- x2aff(X, perplexity, tol = 1e-5, kernel = kernel,
                    verbose = verbose)
    P <- x2ares$W
  }

  # row normalization before anything else
  if (row_normalize) {
    if (symmetrize == "rowsymm") {
      P <- 0.5 * (P + t(P))
      symmetrize <- "none"
    }
    P <- P / rowSums(P)
  }

  # Symmetrize
  P <- switch(symmetrize,
              none = P,
              symmetric = 0.5 * (P + t(P)),
              mutual = sqrt(P * t(P)),
              umap = P + t(P) - P * t(P),
              stop("unknown symmetrization: ", symmetrize))
  # Normalize
  if (normalize) {
    P <- P / sum(P)
  }
  cost$P <- P

  for (r in unique(tolower(ret_extra))) {
    switch(r,
           v = {
             cost$V <- x2ares$W
           },
           dint = {
             if (!is.null(x2ares$dint)) {
               cost$dint <- x2ares$dint
             }
           },
           beta = {
             if (!is.null(x2ares$beta)) {
               cost$beta <- x2ares$beta
             }
           },
           adegc = {
             cost$adegc <- 0.5 * rowSums(x2ares$W) + colSums(x2ares$W)
           },
           adegin = {
             cost$adegin <- rowSums(x2ares$W)
           },
           adegout = {
             cost$adegout <- colSums(x2ares$W)
           },
           pdeg = {
             cost$pdeg <- colSums(P)
           }
    )
  }

  cost
}

# The intrinsic dimensionality associated with a gaussian affinity vector
# Convenient only from in x2aff, where all these values are available
intd_x2aff <- function(D2, beta, W, Z, H, eps = .Machine$double.eps) {
  P <- W / Z
  -2 * beta * sum(D2 * P * (log(P + eps) + H))
}

