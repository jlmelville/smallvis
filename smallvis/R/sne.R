kl_cost <- function(cost, Y) {
  P <- cost$P
  eps <- cost$eps

  cost <- cost_update(cost, Y)

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
      cost <- cost_update(cost, Y)
      P <- cost$P
      cost$G <- k2g(Y, 4 * cost$W * (P - cost$W * cost$invZ))
      cost
    },
    export = function(cost, val) {
      res <- cost_export(cost, val)

      if (is.null(res)) {
        switch(val,
               q = {
                 res <- cost$W * cost$invZ
               })
      }
      res
    },
    update = function(cost, Y) {
      W <- dist2(Y)
      W <- 1 / (1 + W)
      diag(W) <- 0
      invZ <- 1 / sum(W)

      cost$invZ <- invZ
      cost$W <- W
      cost
    },
    sentinel = "W"
  )
}

# Cook, J., Sutskever, I., Mnih, A., & Hinton, G. E. (2007).
# Visualizing similarity data with a mixture of maps.
# In \emph{International Conference on Artificial Intelligence and Statistics} (pp. 67-74).
ssne <- function(perplexity, inp_kernel = "gaussian") {
  lreplace(
    tsne(perplexity = perplexity, inp_kernel = inp_kernel),
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      cost$G <- k2g(Y, 4 * (cost$P - cost$W * cost$invZ))
      cost
    },
    update = function(cost, Y) {
      W <- dist2(Y)
      W <- exp(-W)
      diag(W) <- 0
      invZ <- 1 / sum(W)

      cost$invZ <- invZ
      cost$W <- W
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
      if (is.null(cost$W)) {
        cost <- cost$update(cost, Y)
      }

      invZ <- cost$invZ
      W <- cost$W

      # ASNE defines N KL divergences, row-wise
      cost$pcost <- rowSums(P * log((P + eps) / ((W * invZ) + eps)))
      cost
    },
    gr = function(cost, Y) {
      if (is.null(cost$W)) {
        cost <- cost$update(cost, Y)
      }

      cost$G <- k2g(Y, 2 * (cost$P - cost$W * cost$invZ), symmetrize = TRUE)
      cost
    },
    update = function(cost, Y) {
      W <- dist2(Y)
      W <- exp(-W)
      diag(W) <- 0
      invZ <- 1 / colSums(W)

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
      cost <- cost_update(cost, Y)
      # to include bandwidth
      # K <- 4 * beta * (P - W * invZ) * (W ^ alpha)
      cost$G <- k2g(Y, 4 * (cost$P - cost$W * cost$invZ) * (cost$W ^ alpha))
      cost
    },
    update = function(cost, Y) {
      W <- dist2(Y)
      # to include bandwidth
      # W <- (alpha * beta * W + 1) ^ (-1 / alpha)
      W <- (alpha * W + 1) ^ (-1 / alpha)
      diag(W) <- 0

      cost$invZ <- 1 / sum(W)
      cost$W <- W
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
      cost <- cost_update(cost, Y)
      cost$G <- k2g(Y, 4 * (cost$P - cost$W * cost$invZ) * (cost$W ^ cost$alpha))
      cost
    },
    epoch = function(opt, cost, iter, Y, fn_val) {
      old_cost <- cost

      # create up to two new candidate alpha values:
      # alpha + 0.1 and alpha - 0.1 (but only if that's a valid alpha)
      alpha_curr <- cost$alpha
      alphas <- c(alpha_curr + 0.1)
      if (alpha_curr + 0.1 > alpha_min) {
        alphas <- c(alphas, alpha_curr - 0.1)
      }

      # evaluate error for the new alphas
      errs <- c()
      for (candidate_alpha in alphas) {
        cost$alpha <- candidate_alpha
        cost <- cost_clear(cost)
        cost <- cost_point(cost, Y)
        err <- sum(cost$pcost)
        errs <- c(errs, err)
      }

      # append the current alpha and its error
      errs <- c(errs, fn_val)
      alphas <- c(alphas, alpha_curr)

      # choose the alpha which minimizes the error
      old_cost$alpha <- alphas[which.min(errs)]
      list(cost = old_cost)
    },
    update = function(cost, Y) {
      alpha <- cost$alpha
      W <- dist2(Y)
      W <- (alpha * W + 1) ^ (-1 / alpha)
      diag(W) <- 0

      invZ <- 1 / sum(W)
      cost$invZ <- invZ
      cost$W <- W
      cost
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
      cost <- cost_update(cost, Y)
      cost$G <- k2g(Y, 4 * cost$W * cost$invM * (cost$P - cost$W * cost$invZ))
      cost
    },
    update = function(cost, Y) {
      M <- cost$M

      W <- dist2(Y)
      W <- M / (1 + W)
      diag(W) <- 0
      invZ <- 1 / sum(W)

      cost$invZ <- invZ
      cost$W <- W
      cost
    }
  )
}

wssne <- function(perplexity) {
  lreplace(wtsne(perplexity = perplexity),
     gr = function(cost, Y) {
       cost <- cost_update(cost, Y)
       cost$G <- k2g(Y, 4 * (cost$P - cost$W * cost$invZ))
       cost
     },
     update = function(cost, Y) {
       M <- cost$M

       W <- dist2(Y)
       W <- M * exp(-W)
       diag(W) <- 0
       invZ <- 1 / sum(W)

       cost$invZ <- invZ
       cost$W <- W
       cost
     }
  )
}

# t-SNE but with the gradient defined in terms of un-normalized weights
# Exists entirely as an academic exercise
tsneu <- function(perplexity, inp_kernel = "gaussian") {
  lreplace(
    tsne(perplexity = perplexity, inp_kernel = inp_kernel),
    init = function(cost, X, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      ret_extra = unique(c(ret_extra, "V"))
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = "symmetric", normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra)
      P <- cost$P
      cost$plogp <- colSums(P * log((P + eps)))
      cost$eps <- eps

      cost$V <- cost$V / rowSums(cost$V)
      cost$V <- 0.5 * (cost$V + t(cost$V))
      cost$Vsum <- sum(cost$V)
      cost$invVsum <- 1 / (cost$Vsum)

      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)

      cost$G <- k2g(Y, 4 * cost$W * cost$invVsum * (cost$V - cost$W * cost$invZ * cost$Vsum))
      cost
    }
  )
}

# A pseudo-separable approximation of t-SNE, where the output weight sum is only
# recalculated during the epoch
pstsne <- function(perplexity, inp_kernel = "gaussian") {
  lreplace(
    tsne(perplexity = perplexity, inp_kernel = inp_kernel),
    init = function(cost, X, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      ret_extra = unique(c(ret_extra, "V"))
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = "symmetric", normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra)
      P <- cost$P
      cost$plogp <- colSums(P * log((P + eps)))
      cost$eps <- eps

      # need to row-normalize and symmetrize affinities
      cost$V <- cost$V / rowSums(cost$V)
      cost$V <- 0.5 * (cost$V + t(cost$V))
      cost$Vsum <- sum(cost$V)
      cost$invVsum <- 1 / (cost$Vsum)

      cost$invZ <- 1 / (nrow(cost$P) * nrow(cost$P))

      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)

      cost$G <- k2g(Y, 4 * cost$W * cost$invVsum * (cost$V - cost$W * cost$invZ * cost$Vsum))
      cost
    },
    update = function(cost, Y) {
      W <- dist2(Y)
      W <- 1 / (1 + W)
      diag(W) <- 0

      cost$W <- W

      cost
    },
    epoch = function(opt, cost, iter, Y, fn_val) {
      cost <- cost_update(cost, Y)

      cost$invZ <- 1 / sum(cost$W)
      list(cost = cost)
    }
  )
}

# t-Distributed Elastic Embedding
# EE-like cost function in terms of I-Divergence
# Scaled to give a gradient similar in form to t-SNE
tee <- function(perplexity, inp_kernel = "gaussian", lambda = 0.2) {
  lreplace(
    tsne(perplexity = perplexity, inp_kernel = inp_kernel),
    init = function(cost, X, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      ret_extra = unique(c(ret_extra, "V", "dint"))
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = "symmetric", normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra)
      V <- cost$V
      V <- V / rowSums(V)
      V <- 0.5 * (V + t(V))

      cost$eps <- eps
      cost$V <- V
      cost$invN <- 1 / sum(V)
      cost$gradconst <- 4 * cost$invN
      cost$lambda <- lambda
      cost$constV <- cost$invN * (colSums(V * log(V + eps)) - lambda * colSums(V))

      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      cost$G <- k2g(Y, cost$gradconst * cost$W * (cost$V - cost$W * cost$lambda))
      cost
    },
    update = function(cost, Y) {
      W <- dist2(Y)
      W <- 1 / (1 + W)
      diag(W) <- 0
      cost$W <- W
      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)

      V <- cost$V
      W <- cost$W
      eps <- cost$eps

      cost$pcost <- cost$constV + cost$invN * (cost$lambda * colSums(W) - colSums(V * log(W + eps)))
      cost
    }
  )
}

# UMAP/t-SNE Hybrids ------------------------------------------------------

# Calculate P via normalized smooth knn-distances
skdtsne <- function(perplexity) {
  lreplace(
    tsne(perplexity = perplexity),
    init = function(cost, X, eps = 1e-9, verbose = FALSE, ret_extra = c()) {
      cost$eps <- eps

      P <- smooth_knn_distances(X, k = perplexity, tol = 1e-5,
                                verbose = verbose)$P
      P <- fuzzy_set_union(P)

      # Normalize
      P <- P / sum(P)
      cost$P <- P
      cost$plogp <- colSums(P * log((P + eps)))

      cost
    }
  )
}

# Use the UMAP curve family in output kernel
usne <- function(perplexity, inp_kernel = "gaussian", spread = 1,
                  min_dist = 0.001, gr_eps = 0.1) {
  lreplace(
    tsne(perplexity = perplexity),
    init = function(cost, X, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = "symmetric", normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra)
      P <- cost$P
      cost$plogp <- colSums(P * log((P + eps)))
      cost$eps <- eps

      ab_params <- find_ab_params(spread = spread, min_dist = min_dist)
      a <- ab_params[1]
      b <- ab_params[2]
      if (verbose) {
        message("Umap curve parameters = ", formatC(a), ", ", formatC(b))
      }
      cost$a <- a
      cost$b <- b

      cost
    },
    update = function(cost, Y) {
      D2 <- dist2(Y)
      D2[D2 < 0] <- 0

      W <- 1 / (1 + cost$a * D2 ^ cost$b)
      diag(W) <- 0

      cost$invZ <- 1 / sum(W)
      cost$W <- W
      cost$D2 <- D2
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      cost$G <- k2g(Y, 4 * (cost$b * (1 - cost$W) / (cost$D2 + gr_eps)) * (cost$P - cost$W * cost$invZ))
      cost
    }
  )
}

# UMAP cross entropy cost instead of KL divergence
cetsne <- function(perplexity, inp_kernel = "gaussian") {
  lreplace(tsne(perplexity),
           init = function(cost, X, eps = .Machine$double.eps, verbose = FALSE,
                           ret_extra = c()) {
             cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                              symmetrize = "symmetric", normalize = TRUE,
                              verbose = verbose, ret_extra = ret_extra)
             P <- cost$P
             cost$Cp <- colSums(P * log(P + eps) + (1 - P) * log1p(-P + eps))
             cost$eps <- eps
             cost
           },
           pfn = function(cost, Y) {
             cost <- cost_update(cost, Y)

             P <- cost$P
             eps <- cost$eps
             Q <- cost$Q

             cost$pcost <- colSums(-P * log(Q + eps) - (1 - P) * log1p(-Q + eps)) + cost$Cp
             cost
           },
           gr = function(cost, Y) {
             cost <- cost_update(cost, Y)
             cost$G <- k2g(Y,  4 * cost$W * (cost$C - cost$sumC * cost$Q))
             cost
           },
           update = function(cost, Y) {
             P <- cost$P
             W <- dist2(Y)
             W <- 1 / (1 + W)
             diag(W) <- 0

             invZ <- 1 / sum(W)
             Q <- W * invZ
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
      cost <- cost_update(cost, Y)
      cost$G <- k2g(Y, 2 * cost$beta * cost$W * (cost$P - cost$W * cost$invZ), symmetrize = TRUE)
      cost
    },
    update = function(cost, Y) {
      W <- dist2(Y)
      W <- 1 / (1 + (cost$beta * W))
      diag(W) <- 0

      cost$invZ <- 1 / sum(W)
      cost$W <- W
      cost
    }
  )
}

# SSNE with input kernel bandwidths transferred to output
bssne <- function(perplexity, inp_kernel = "gaussian") {
  lreplace(btsne(perplexity = perplexity, inp_kernel = inp_kernel),
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      cost$G <- k2g(Y, 2 * cost$beta * (cost$P - cost$W * cost$invZ), symmetrize = TRUE)
      cost
    },
    update = function(cost, Y) {
      W <- dist2(Y)
      W <- exp(-cost$beta * W)
      diag(W) <- 0
      invZ <- 1 / sum(W)

      cost$invZ <- 1 / sum(W)
      cost$W <- W
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
    cost <- cost_update(cost, Y)

    cost$G <- k2g(Y, 2 * cost$W * (cost$P - cost$W * cost$invZ), symmetrize = TRUE)
    cost
  },
  update = function(cost, Y) {
    W <- dist2(Y)
    W <- 1 / (1 + W)
    diag(W) <- 0

    cost$W <- W
    cost$invZ <- 1 / rowSums(W)
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
    if (is.character(perplexity) || is.list(perplexity)) {
      stop("Can't use intrinsic dimensionality with knn kernel")
    }
    if (length(perplexity) > 1) {
      stop("Can't use multiple perplexities with knn kernel")
    }
    if (verbose) {
      tsmessage("Using knn kernel with k = ", formatC(perplexity))
    }
    P <- knn_graph(X, k = perplexity)
    x2ares <- list(W = P)
  }
  if (perp_method(perplexity) == "idp") {
    perplexities <- NULL
    if (is.list(perplexity) && length(perplexity) == 2) {
      perplexities <- perplexity[[2]]
    }

    x2ares <- idp(X, perplexities = perplexities, tol = 1e-5,
                  verbose = verbose)
    P <- x2ares$W
    ret_extra <- unique(c(ret_extra, "idp"))
  }
  else {
    if (!is.numeric(perplexity)) {
      stop("Unknown perplexity method, '", perplexity[[1]], "'")
    }
    if (verbose) {
      tsmessage("Commencing calibration for perplexity = ",
                format_perps(perplexity))
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
           },
           idp = {
             if (!is.null(x2ares$idp)) {
              cost$idp <- x2ares$idp
             }
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

