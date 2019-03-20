# Calculates shifted exponential column-wise: exp(X - a)
# where a is the column max.
# This is the log-sum-exp trick to avoid numeric underflow:
# log sum_i exp x_i = a + log sum_i exp(x_i - a)
# => sum_i exp x_i = exp a * sum_i exp(x_i - a)
# with a = max x_i
# exp(max x_i) can still underflow so we don't return Z (the sum)
# Use Q directly (exp a appears in numerator and denominator, so cancels).
# https://statmodeling.stat.columbia.edu/2016/06/11/log-sum-of-exponentials/
# https://www.xarg.org/2016/06/the-log-sum-exp-trick-in-machine-learning/
# http://wittawat.com/posts/log-sum_exp_underflow.html
exp_shift <- function(X) {
  X <- exp(sweep(X, 2, apply(X, 2, max)))
}

expQ <- function(Y, eps = .Machine$double.xmin, beta = NULL,
                     A = NULL,
                     is_symmetric = FALSE,
                     matrix_normalize = FALSE) {
  W <- dist2(Y)
  
  if (!is.null(beta)) {
    W <- exp_shift(-W * beta)
  }
  else {
    W <- exp_shift(-W)
  }
  
  if (!is.null(A)) {
    W <- A * W
  }
  diag(W) <- 0
  
  if (matrix_normalize) {
    Z <- sum(W)
  }
  else {
    if (is_symmetric) {
      Z <- colSums(W)
    }
    else {
      Z <- rowSums(W)
    }
  }
  # cost of division (vs storing 1/Z and multiplying) seems small
  Q <- W / Z
  
  if (eps > 0) {
    Q[Q < eps] <- eps
  }
  diag(Q) <- 0
  
  list(
    Q = Q,
    Z = Z
  )
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

# t-SNE
tsne <- function(perplexity, inp_kernel = "gaussian") {
  list(
    init = function(cost, X, max_iter, eps = .Machine$double.xmin, verbose = FALSE,
                    ret_extra = c()) {
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = "symmetric", normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra)
      cost$eps <- eps
      cost
    },
    cache_input = function(cost) {
      # cache P log P constant part of cost: incur 1 extra log operation now
      # but saves one division operation every time we calculate cost:
      # substantial (10-15%) speed up with Wolfe line search methods
      P <- cost$P
      eps <- cost$eps
      P[P < eps] <- eps
      cost$plogp <- colSums(P * logm(P, eps))
      cost
    },
    pfn = kl_cost,
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      P <- cost$P
      cost$G <- k2g(Y, 4 * cost$W * (P - cost$W / cost$Z))
      
      cost
    },
    export = function(cost, val) {
      res <- cost_export(cost, val)

      if (is.null(res)) {
        switch(val,
               q = {
                 res <- cost$W / cost$Z
               })
      }
      res
    },
    update = function(cost, Y) {
      W <- dist2(Y)
      W <- 1 / (1 + W)
      diag(W) <- 0

      cost$Z <- sum(W)
      cost$W <- W
      cost
    },
    sentinel = "W",
    exaggerate = function(cost, exaggeration_factor) {
      cost$P <- cost$P * exaggeration_factor
      cost
    }
  )
}

# Cook, J., Sutskever, I., Mnih, A., & Hinton, G. E. (2007).
# Visualizing similarity data with a mixture of maps.
# In \emph{International Conference on Artificial Intelligence and Statistics} (pp. 67-74).
ssne <- function(perplexity, inp_kernel = "gaussian") {
  lreplace(
    tsne(perplexity = perplexity, inp_kernel = inp_kernel),
    pfn = kl_costQ,
    gr = function(cost, Y) {
      cost <- cost$update(cost, Y)
      cost$G <- k2g(Y, 4 * (cost$P - cost$Q), symmetrize = FALSE)
      
      cost
    },
    update = function(cost, Y) {
      cost$Q <- expQ(Y, cost$eps, is_symmetric = TRUE,
                         matrix_normalize = TRUE)$Q
      cost
    },
    sentinel = "Q"
  )
}

# Hinton, G. E., & Roweis, S. T. (2002).
# Stochastic neighbor embedding.
# In \emph{Advances in neural information processing systems} (pp. 833-840).
asne <- function(perplexity, inp_kernel = "gaussian") {
  lreplace(tsne(perplexity),
    init = function(cost, X, max_iter, eps = .Machine$double.xmin, verbose = FALSE,
                    ret_extra = c()) {
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = "none", normalize = FALSE,
                       verbose = verbose, ret_extra = ret_extra)
      cost$eps <- eps
      cost
    },
    pfn = kl_costQr,
    gr = function(cost, Y) {
      cost <- cost$update(cost, Y)
      cost$G <- k2g(Y, 2 * (cost$P - cost$Q), symmetrize = TRUE)

      cost
    },
    update = function(cost, Y) {
      cost$Q <- expQ(Y, eps = cost$eps, is_symmetric = FALSE)$Q
      cost
    },
    sentinel = "Q"
  )
}

# Heavy-Tailed Symmetric Stochastic Neighbor Embedding (HSSNE)
# Yang, Z., King, I., Xu, Z., & Oja, E. (2009).
# Heavy-tailed symmetric stochastic neighbor embedding.
# In \emph{Advances in neural information processing systems} (pp. 2169-2177).
hssne <- function(perplexity, alpha = 0.5) {
  alpha <- max(alpha, 1e-8)
  apow <- -1 / alpha
  lreplace(
    tsne(perplexity = perplexity),
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      # to include bandwidth
      # K <- 4 * beta * (P - W / Z) * powm(W, alpha, eps)
      W <- cost$W
      cost$G <- k2g(Y, 4 * (cost$P - W / cost$Z) * powm(W, alpha, cost$eps))
      cost
    },
    update = function(cost, Y) {
      W <- dist2(Y)
      # to include bandwidth
      # W <- (alpha * beta * W + 1) ^ (-1 / alpha)
      W <- powm(alpha * W + 1, apow, cost$eps)
      diag(W) <- 0

      cost$Z <- sum(W)
      cost$W <- W
      cost
    }
  )
}

# exists to demonstrate that contant beta doesn't have any meaningful effect
# on the results.
bhssne <- function(perplexity, alpha = 0.5, beta = 1) {
  alpha <- max(alpha, 1e-8)
  beta <- max(beta, 1e-8)
  lreplace(
    tsne(perplexity = perplexity),
    b4 = 4 * beta,
    ab = alpha * beta,
    apow = -1 / alpha,
    init = function(cost, X, max_iter, eps = .Machine$double.xmin, verbose = FALSE,
                    ret_extra = c()) {
      ret_extra <- unique(c(ret_extra, 'beta'))
      cost <- sne_init(cost, X, perplexity = perplexity,
                       symmetrize = "symmetric", normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra)
      # override input bandwidths with fixed beta (although this doesn't do much)
      if (!is.null(beta)) {
        cost$beta <- beta
      }
      cost$eps <- eps
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      W <- cost$W
      cost$G <- k2g(Y, cost$b4 * (cost$P - W / cost$Z) * powm(W, alpha, cost$eps))

      cost
    },
    update = function(cost, Y) {
      W <- dist2(Y)
      W <- powm(cost$ab * W + 1, cost$apow, cost$eps)
      diag(W) <- 0
      
      cost$Z <- sum(W)
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
      cost$G <- k2g(Y, 4 * (cost$P - cost$W / cost$Z) * powm(cost$W, cost$alpha, cost$eps))
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
        cost$apow <- -1 / candidate_alpha
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
      old_cost$apow <- -1 / old_cost$alpha
      list(cost = old_cost)
    },
    update = function(cost, Y) {
      alpha <- cost$alpha
      W <- dist2(Y)
      W <- powm(alpha * W + 1, cost$apow, cost$eps)
      diag(W) <- 0

      cost$Z <- sum(W)
      cost$W <- W
      cost
    },
    alpha = alpha,
    apow = -1 / alpha
  )
}

# Yang, Z., Peltonen, J., & Kaski, S. (2014).
# Optimization equivalence of divergences improves neighbor embedding.
# In \emph{Proceedings of the 31st International Conference on Machine Learning (ICML-14)}
# (pp. 460-468).
wtsne <- function(perplexity) {
  lreplace(tsne(perplexity = perplexity),
    init = function(cost, X, max_iter, eps = .Machine$double.xmin, verbose = FALSE,
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

      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      cost$G <- k2g(Y, 4 * cost$W * cost$invM * (cost$P - cost$W / cost$Z))
      cost
    },
    update = function(cost, Y) {
      M <- cost$M

      W <- dist2(Y)
      W <- M / (1 + W)
      diag(W) <- 0

      cost$Z <- sum(W)
      cost$W <- W
      cost
    }
  )
}

wssne <- function(perplexity) {
  lreplace(ssne(perplexity = perplexity),
     init = function(cost, X, max_iter, eps = .Machine$double.xmin, verbose = FALSE,
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
       
       cost
     },
     gr = function(cost, Y) {
       cost <- cost_update(cost, Y)
       cost$G <- k2g(Y, 4 * (cost$P - cost$Q))
       cost
     },
     update = function(cost, Y) {
       cost$Q <- expQ(Y, cost$eps, A = cost$M, matrix_normalize = TRUE)$Q
       cost
     }
  )
}

# t-SNE but with the gradient defined in terms of un-normalized weights
# Exists entirely as an academic exercise
tsneu <- function(perplexity, inp_kernel = "gaussian") {
  lreplace(
    tsne(perplexity = perplexity, inp_kernel = inp_kernel),
    init = function(cost, X, max_iter, eps = .Machine$double.xmin, verbose = FALSE,
                    ret_extra = c()) {
      ret_extra = unique(c(ret_extra, "V"))
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = "symmetric", normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra)

      cost$eps <- eps

      cost$V <- cost$V / rowSums(cost$V)
      cost$V <- 0.5 * (cost$V + t(cost$V))
      cost$Vsum <- sum(cost$V)
      cost$invVsum <- 1 / (cost$Vsum)

      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)

      cost$G <- k2g(Y, 4 * cost$W * cost$invVsum * (cost$V - (cost$W / cost$Z) * cost$Vsum))
      cost
    }
  )
}

# A pseudo-separable approximation of t-SNE, where the output weight sum is only
# recalculated during the epoch
pstsne <- function(perplexity, inp_kernel = "gaussian") {
  lreplace(
    tsne(perplexity = perplexity, inp_kernel = inp_kernel),
    init = function(cost, X, max_iter, eps = .Machine$double.xmin, verbose = FALSE,
                    ret_extra = c()) {
      ret_extra = unique(c(ret_extra, "V"))
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = "symmetric", normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra)

      # need to row-normalize and symmetrize affinities
      cost$V <- cost$V / rowSums(cost$V)
      cost$V <- 0.5 * (cost$V + t(cost$V))
      cost$Vsum <- sum(cost$V)
      cost$invVsum <- 1 / (cost$Vsum)

      cost
    },
    cache_input = function(cost) {
      P <- cost$P
      eps <- cost$eps
      P[P < eps] <- eps
      cost$plogp <- colSums(P * logm(P, eps))
      cost$Z <- nrow(P) * nrow(P)
      
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)

      cost$G <- k2g(Y, 4 * cost$W * cost$invVsum * (cost$V - (cost$W / cost$Z) * cost$Vsum))
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

      cost$Z <- sum(cost$W)
      list(cost = cost)
    }
  )
}

# t-Distributed Elastic Embedding
# EE-like cost function in terms of I-Divergence
# Scaled to give a gradient similar in form to t-SNE
tee <- function(perplexity, inp_kernel = "gaussian", lambda = 0.01) {
  lreplace(
    tsne(perplexity = perplexity, inp_kernel = inp_kernel),
    init = function(cost, X, max_iter, eps = .Machine$double.xmin, verbose = FALSE,
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
      
      V[V < eps] <- eps
      cost$constV <- cost$invN * (colSums(V * logm(V, eps)) - lambda * colSums(V))
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

      cost$pcost <- cost$constV +
        cost$invN * (cost$lambda * colSums(W) - colSums(V * logm(W, eps)))
      cost
    }
  )
}

# UMAP/t-SNE Hybrids ------------------------------------------------------

# Calculate P via normalized smooth knn-distances
skdtsne <- function(perplexity) {
  lreplace(
    tsne(perplexity = perplexity),
    init = function(cost, X, max_iter, eps = 1e-9, verbose = FALSE, ret_extra = c()) {
      cost$eps <- eps

      P <- smooth_knn_distances(X, k = perplexity, tol = 1e-5,
                                verbose = verbose)$P
      P <- fuzzy_set_union(P)
      # Normalize
      P <- P / sum(P)
      cost$P <- P

      cost
    }
  )
}

# Use the UMAP curve family in output kernel
usne <- function(perplexity, inp_kernel = "gaussian", spread = 1,
                  min_dist = 0.001, gr_eps = 0.1) {
  lreplace(
    tsne(perplexity = perplexity),
    init = function(cost, X, max_iter, eps = .Machine$double.xmin, verbose = FALSE,
                    ret_extra = c()) {
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = "symmetric", normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra)
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

      cost$Z <- sum(W)
      cost$W <- W
      cost$D2 <- D2
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      cost$G <- k2g(Y, 4 * (cost$b * (1 - cost$W) / (cost$D2 + gr_eps)) * (cost$P - cost$W / cost$Z))
      cost
    }
  )
}

# UMAP cross entropy cost instead of KL divergence
cetsne <- function(perplexity, inp_kernel = "gaussian") {
  lreplace(tsne(perplexity),
           init = function(cost, X, max_iter, eps = .Machine$double.xmin, verbose = FALSE,
                           ret_extra = c()) {
             cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                              symmetrize = "symmetric", normalize = TRUE,
                              verbose = verbose, ret_extra = ret_extra)

             cost$eps <- eps
             cost
           },
           cache_input = function(cost) {
             P <- cost$P
             eps <- cost$eps
             P[P < eps] <- eps
             cost$Cp <- colSums(P * logm(P, eps) + (1 - P) * log1p(-P))
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
             cost$G <- k2g(Y,  4 * cost$W * (cost$C - cost$sumC * cost$Q))
             cost
           },
           update = function(cost, Y) {
             P <- cost$P
             W <- dist2(Y)
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


# Bandwidth Experiments --------------------------------------------------

# t-SNE with input kernel bandwidths transferred to output
btsne <- function(perplexity, inp_kernel = "gaussian", beta = NULL) {
  lreplace(tsne(perplexity = perplexity, inp_kernel = inp_kernel),
    init = function(cost, X, max_iter, eps = .Machine$double.xmin, verbose = FALSE,
                    ret_extra = c()) {
      ret_extra <- unique(c(ret_extra, 'beta'))
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = "symmetric", normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra)
      # override input bandwidths with fixed beta (although this doesn't do much)
      if (!is.null(beta)) {
        cost$beta <- beta
      }
      cost$eps <- eps
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      cost$G <- k2g(Y, 2 * cost$beta * cost$W * (cost$P - (cost$W / cost$Z)), symmetrize = TRUE)
      cost
    },
    update = function(cost, Y) {
      W <- dist2(Y)
      W <- 1 / (1 + (cost$beta * W))
      diag(W) <- 0

      cost$Z <- sum(W)
      cost$W <- W
      cost
    }
  )
}

# SSNE with input kernel bandwidths transferred to output
bssne <- function(perplexity, inp_kernel = "gaussian", beta = NULL) {
  lreplace(ssne(perplexity = perplexity, inp_kernel = inp_kernel),
    init = function(cost, X, max_iter, eps = .Machine$double.xmin, verbose = FALSE,
                     ret_extra = c()) {
       ret_extra <- unique(c(ret_extra, 'beta'))
       cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                        symmetrize = "symmetric", normalize = TRUE,
                        verbose = verbose, ret_extra = ret_extra)
       if (!is.null(beta)) {
         cost$beta <- beta
       }
       cost$eps <- eps
       cost
    },
    pfn = kl_costQ,
    gr = function(cost, Y) {
      cost <- cost$update(cost, Y)
      cost$G <- k2g(Y, 2 * cost$beta * (cost$P - cost$Q), symmetrize = TRUE)
      
      cost
    },
    update = function(cost, Y) {
      cost$Q <- expQ(Y, cost$eps, beta = cost$beta, matrix_normalize = TRUE)$Q
      cost
    }
  )
}

# ASNE with input kernel bandwidths transferred to output
basne <- function(perplexity, beta = NULL) {
  lreplace(
    asne(perplexity = perplexity),
    init = function(cost, X, max_iter, eps = .Machine$double.xmin, verbose = FALSE,
                    ret_extra = c()) {
      ret_extra <- unique(c(ret_extra, 'beta'))

      cost <- sne_init(cost, X, perplexity = perplexity,
                       symmetrize = "none", normalize = FALSE,
                       verbose = verbose, ret_extra = ret_extra)
      # override input bandwidths with fixed beta (although this doesn't do much)
      if (!is.null(beta)) {
        cost$beta <- beta
      }

      cost$eps <- eps
      cost
    },
    update = function(cost, Y) {
      cost$Q <- expQ(Y, eps = cost$eps, beta = cost$beta, is_symmetric = FALSE)$Q
      cost
    },
    gr = function(cost, Y) {
      cost <- cost$update(cost, Y)
      cost$G <- k2g(Y, 2 * cost$beta * (cost$P - cost$Q), symmetrize = TRUE)
      cost
   }
  )
}

# t-ASNE with input kernel bandwidths transferred to output
btasne <- function(perplexity, beta = NULL) {
  lreplace(basne(
    perplexity = perplexity, beta = beta),
    pfn = kl_cost,
    gr = function(cost, Y) {
    cost <- cost_update(cost, Y)
     W <- cost$W
     cost$G <- k2g(Y, 2 * cost$beta * W * (cost$P - W / cost$Z), 
                   symmetrize = TRUE)
     cost
    },
    update = function(cost, Y) {
      W <- dist2(Y)
      W <- 1 / (1 + cost$beta * W)
      diag(W) <- 0
       
      cost$W <- W
      cost$Z <- rowSums(W)
      cost
    }
  )
}

tasne <- function(perplexity) {
  lreplace(tsne(perplexity = perplexity),
           init = function(cost, X, max_iter, eps = .Machine$double.xmin, verbose = FALSE,
                           ret_extra = c()) {
             cost <- sne_init(cost, X, perplexity = perplexity,
                              symmetrize = "none", normalize = FALSE,
                              verbose = verbose, ret_extra = ret_extra)
             cost$eps <- eps
             cost
           },
           gr = function(cost, Y) {
             cost <- cost_update(cost, Y)
             
             cost$G <- k2g(Y, 2 * cost$W * (cost$P - cost$W / cost$Z), 
                           symmetrize = TRUE)
             cost
           },
           update = function(cost, Y) {
             W <- dist2(Y)
             W <- 1 / (1 + W)
             diag(W) <- 0
             
             cost$W <- W
             cost$Z <- 1 / rowSums(W)
             cost
           }
  )
}


# Normalization Experiments -----------------------------------------------

# ASNE but with the t-distributed kernel
tasne <- function(perplexity) {
  lreplace(tsne(perplexity = perplexity),
  init = function(cost, X, max_iter, eps = .Machine$double.xmin, verbose = FALSE,
                  ret_extra = c()) {
    cost <- sne_init(cost, X, perplexity = perplexity,
                     symmetrize = "none", normalize = FALSE,
                     verbose = verbose, ret_extra = ret_extra)
    cost$eps <- eps
    cost
  },
  gr = function(cost, Y) {
    cost <- cost_update(cost, Y)

    cost$G <- k2g(Y, 2 * cost$W * (cost$P - cost$W / cost$Z), 
                  symmetrize = TRUE)
    cost
  },
  update = function(cost, Y) {
    W <- dist2(Y)
    W <- 1 / (1 + W)
    diag(W) <- 0

    cost$W <- W
    cost$Z <- rowSums(W)
    cost
  }
  )
}

# t-RM-SNE
# t-SNE without symmetrization of P (but still pair-normalizing)
# row-normalize, then matrix normalize
trmsne <- function(perplexity, inp_kernel = "gaussian") {
  lreplace(tsne(perplexity = perplexity, inp_kernel = inp_kernel),
           init = function(cost, X, max_iter, eps = .Machine$double.xmin, verbose = FALSE,
                           ret_extra = c()) {
             cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                              symmetrize = "none", normalize = TRUE,
                              verbose = verbose, ret_extra = ret_extra)
             cost$eps <- eps
             cost
           },
           gr = function(cost, Y) {
             cost <- cost_update(cost, Y)

             W <- cost$W
             cost$G <- k2g(Y, 2 * W * (cost$P - W / cost$Z), symmetrize = TRUE)
             cost
           }
  )
}

# t-M-SNE
# t-SNE but without row-normalizing or symmetrizing, just matrix normalization
# Not recommended
tmsne <- function(perplexity, inp_kernel = "gaussian") {
  lreplace(trmsne(perplexity = perplexity, inp_kernel = inp_kernel),
           init = function(cost, X, max_iter, eps = .Machine$double.xmin, verbose = FALSE,
                           ret_extra = c()) {
             cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                              symmetrize = "none", row_normalize = FALSE,
                              normalize = TRUE,
                              verbose = verbose, ret_extra = ret_extra)
             cost$eps <- eps
             cost
           }
  )
}

# RSR row-normalize, symmetrize, then row-normalize again
# Might work a tiny bit better than t-ASNE?
trsrsne <- function(perplexity) {
  lreplace(tasne(perplexity),
           init = function(cost, X, max_iter, eps = .Machine$double.xmin, verbose = FALSE,
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
           init = function(cost, X, max_iter, eps = .Machine$double.xmin, verbose = FALSE,
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

scale_affinities <- function(P, symmetrize = "symmetric", row_normalize = TRUE,
                             normalize = TRUE) {
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
  P
}


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
  else if (perp_method(perplexity) == "multiscale") {
    perplexities <- NULL
    if (is.list(perplexity) && length(perplexity) == 2) {
      perplexities <- perplexity[[2]]
    }
    
    mspres <- msp(X, perplexities = perplexities, tol = 1e-5,
                  symmetrize = symmetrize, 
                  row_normalize = row_normalize,
                  normalize = normalize,
                  verbose = verbose)
    cost$P <- mspres$P
    return(cost)
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

  P <- scale_affinities(P, 
                        symmetrize = symmetrize, 
                        row_normalize = row_normalize,
                        normalize = normalize)
  
  cost$P <- P

  tsmessage("Effective perplexity of P approx = ", 
            formatC(stats::median(perpp(P))))
  
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
intd_x2aff <- function(D2, beta, W, Z, H, eps = .Machine$double.xmin) {
  P <- W / Z
  -2 * beta * sum(D2 * P * (log(P + eps) + H))
}

shannonpr <- function(P, eps = .Machine$double.xmin) {
  P <- P / rowSums(P)
  rowSums(-P * log(P + eps))
}

perpp <- function(P) {
  exp(shannonpr(P))
}

