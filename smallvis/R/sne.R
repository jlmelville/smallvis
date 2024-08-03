# t-SNE
tsne <- function(perplexity, inp_kernel = "gaussian", symmetrize = "symmetric", 
                 normalize = TRUE, row_normalize = TRUE,
                 eps = .Machine$double.eps, n_threads = 0, use_cpp = FALSE) {
  list(
    init = function(cost, X, max_iter, verbose = FALSE,
                    ret_extra = c()) {
      symmetrize <- match.arg(tolower(symmetrize), true_symmetrize_options())
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = symmetrize, normalize = normalize,
                       row_normalize = row_normalize,
                       verbose = verbose, ret_extra = ret_extra,
                       n_threads = n_threads, use_cpp = use_cpp)
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
      if (use_cpp) {
        cost$G <- tsne_grad_cpp(P, cost$W, cost$Z, Y, n_threads = n_threads)
      }
      else {
        cost$G <- k2g(Y, 4 * cost$W * (P - cost$W / cost$Z))
      }
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
      W <- calc_tweight(Y, use_cpp = use_cpp, n_threads = n_threads)
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
ssne <- function(perplexity, inp_kernel = "gaussian", 
                 symmetrize = "symmetric", eps = .Machine$double.eps,
                 n_threads = 0, use_cpp = FALSE) {
  lreplace(
    tsne(perplexity = perplexity, inp_kernel = inp_kernel, 
         symmetrize = symmetrize, eps = eps, n_threads = n_threads,
         use_cpp = use_cpp),
    pfn = kl_costQ,
    gr = function(cost, Y) {
      cost <- cost$update(cost, Y)
      cost$G <- k2g(Y, 4 * (cost$P - cost$Q), symmetrize = FALSE)
      
      cost
    },
    update = function(cost, Y) {
      cost$Q <- expQ(Y, cost$eps, is_symmetric = TRUE, matrix_normalize = TRUE,
                     use_cpp = use_cpp, n_threads = n_threads)$Q
      cost
    },
    sentinel = "Q"
  )
}

# Hinton, G. E., & Roweis, S. T. (2002).
# Stochastic neighbor embedding.
# In \emph{Advances in neural information processing systems} (pp. 833-840).
asne <- function(perplexity, inp_kernel = "gaussian", 
                 eps = .Machine$double.eps, n_threads = 0, use_cpp = FALSE) {
  lreplace(tsne(perplexity, use_cpp = use_cpp, n_threads = n_threads),
    init = function(cost, X, max_iter, verbose = FALSE,
                    ret_extra = c()) {
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = "none", normalize = FALSE,
                       verbose = verbose, ret_extra = ret_extra,
                       n_threads = n_threads, use_cpp = use_cpp)
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
      cost$Q <- expQ(Y, eps = cost$eps, is_symmetric = FALSE,
                     use_cpp = use_cpp, n_threads = n_threads)$Q
      cost
    },
    sentinel = "Q"
  )
}

# Heavy-Tailed Symmetric Stochastic Neighbor Embedding (HSSNE)
# Yang, Z., King, I., Xu, Z., & Oja, E. (2009).
# Heavy-tailed symmetric stochastic neighbor embedding.
# In \emph{Advances in neural information processing systems} (pp. 2169-2177).
hssne <- function(perplexity, alpha = 0.5, inp_kernel = "gaussian", 
                  symmetrize = "symmetric", eps = .Machine$double.eps,
                  n_threads = 0, use_cpp = FALSE) {
  alpha <- max(alpha, 1e-8)
  apow <- -1 / alpha
  lreplace(
    tsne(perplexity = perplexity, inp_kernel = inp_kernel, 
         symmetrize = symmetrize, eps = eps, n_threads = n_threads),
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      # to include bandwidth
      # K <- 4 * beta * (P - W / Z) * powm(W, alpha, eps)
      W <- cost$W
      cost$G <- k2g(Y, 4 * (cost$P - W / cost$Z) * powm(W, alpha, cost$eps))
      cost
    },
    update = function(cost, Y) {
      W <- calc_d2(Y, use_cpp = use_cpp, n_threads = n_threads)
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

# exists to demonstrate that constant beta doesn't have any meaningful effect
# on the results.
bhssne <- function(perplexity, alpha = 0.5, beta = 1,
                   eps = .Machine$double.eps, n_threads = 0, use_cpp = FALSE) {
  alpha <- max(alpha, 1e-8)
  beta <- max(beta, 1e-8)
  lreplace(
    tsne(perplexity = perplexity, use_cpp = use_cpp, n_threads = n_threads),
    b4 = 4 * beta,
    ab = alpha * beta,
    apow = -1 / alpha,
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      ret_extra <- unique(c(ret_extra, 'beta'))
      cost <- sne_init(cost, X, perplexity = perplexity,
                       symmetrize = "symmetric", normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra,
                       n_threads = n_threads, use_cpp = use_cpp)
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
      W <- calc_d2(Y, use_cpp = use_cpp, n_threads = n_threads)
      W <- powm(cost$ab * W + 1, cost$apow, cost$eps)
      diag(W) <- 0
      
      cost$Z <- sum(W)
      cost$W <- W
      cost
    }
  )
}

# A version of HSSNE where alpha is allowed to vary at every epoch
dhssne <- function(perplexity, alpha = 0.5, inp_kernel = "gaussian", 
                   symmetrize = "symmetric", eps = .Machine$double.eps,
                   n_threads = 0, use_cpp = FALSE) {
  alpha_min <- 1e-8
  alpha <- max(alpha, alpha_min)
  lreplace(
    tsne(perplexity = perplexity, inp_kernel = inp_kernel, 
         symmetrize = symmetrize, eps = eps, n_threads = n_threads),
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
      W <- calc_d2(Y, use_cpp = use_cpp, n_threads = n_threads)
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
wtsne <- function(perplexity, inp_kernel = "gaussian", 
                  symmetrize = "symmetric", eps = .Machine$double.eps,
                  n_threads = 0, use_cpp = FALSE) {
  lreplace(tsne(perplexity = perplexity, use_cpp = use_cpp,
                n_threads = n_threads),
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      symmetrize <- match.arg(tolower(symmetrize), true_symmetrize_options())
      ret_extra <- c(ret_extra, "pdeg")
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = symmetrize, normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra,
                       n_threads = n_threads, use_cpp = use_cpp)
      # P matrix degree centrality: column sums
      deg <- cost$pdeg
      if (verbose) {
        summarize(deg, "deg", verbose = verbose)
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

      W <- calc_d2(Y, use_cpp = use_cpp, n_threads = n_threads)
      W <- M / (1 + W)
      diag(W) <- 0

      cost$Z <- sum(W)
      cost$W <- W
      cost
    }
  )
}

wssne <- function(perplexity, inp_kernel = "gaussian", 
                  symmetrize = "symmetric", eps = .Machine$double.eps,
                  n_threads = 0, use_cpp = FALSE) {
  lreplace(ssne(perplexity = perplexity, use_cpp = use_cpp,
                n_threads = n_threads),
     init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
       symmetrize <- match.arg(tolower(symmetrize), true_symmetrize_options())
       ret_extra <- c(ret_extra, "pdeg")
       cost <- sne_init(cost, X, perplexity = perplexity,
                        kernel = inp_kernel, symmetrize = symmetrize, 
                        normalize = TRUE, verbose = verbose, 
                        ret_extra = ret_extra, n_threads = n_threads,
                        use_cpp = use_cpp)
       # P matrix degree centrality: column sums
       deg <- cost$pdeg
       if (verbose) {
         summarize(deg, "deg", verbose = verbose)
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
       cost$Q <- expQ(Y, cost$eps, A = cost$M, matrix_normalize = TRUE, 
                      use_cpp = use_cpp, n_threads = n_threads)$Q
       cost
     }
  )
}

# t-SNE but with the gradient defined in terms of un-normalized weights
# Exists entirely as an academic exercise
tsneu <- function(perplexity, inp_kernel = "gaussian", 
                  eps = .Machine$double.eps, n_threads = 0, use_cpp = FALSE) {
  lreplace(
    tsne(perplexity = perplexity, inp_kernel = inp_kernel, use_cpp = use_cpp,
         n_threads = n_threads),
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      ret_extra = unique(c(ret_extra, "V"))
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = "symmetric", normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra,
                       n_threads = n_threads, use_cpp = use_cpp)

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
pstsne <- function(perplexity, inp_kernel = "gaussian",
                   eps = .Machine$double.eps, n_threads = 0, use_cpp = FALSE) {
  lreplace(
    tsne(perplexity = perplexity, inp_kernel = inp_kernel, use_cpp = use_cpp,
         n_threads = n_threads),
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      ret_extra = unique(c(ret_extra, "V"))
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = "symmetric", normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra,
                       n_threads = n_threads, use_cpp = use_cpp)

      # need to row-normalize and symmetrize affinities
      cost$V <- cost$V / rowSums(cost$V)
      cost$V <- 0.5 * (cost$V + t(cost$V))
      cost$Vsum <- sum(cost$V)
      cost$invVsum <- 1 / (cost$Vsum)

      cost$eps <- eps
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
      W <- calc_d2(Y, use_cpp = use_cpp, n_threads = n_threads)
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

# UMAP/t-SNE Hybrids ------------------------------------------------------

# Calculate P via normalized smooth knn-distances
skdtsne <- function(perplexity, eps = .Machine$double.eps, n_threads = 0,
                    use_cpp = FALSE) {
  tsne(perplexity = perplexity, inp_kernel = "skd", symmetrize = "umap", 
       eps = eps, n_threads = n_threads, use_cpp = use_cpp)
}

# Use the UMAP curve family in output kernel
usne <- function(perplexity, inp_kernel = "gaussian", symmetrize = "symmetric", 
                 spread = 1, min_dist = 0.001, gr_eps = 0.1, 
                 eps = .Machine$double.eps, n_threads = 0, use_cpp = FALSE) {
  lreplace(
    tsne(perplexity = perplexity, use_cpp = use_cpp, n_threads = n_threads),
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      symmetrize <- match.arg(tolower(symmetrize), true_symmetrize_options())
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = symmetrize, normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra,
                       n_threads = n_threads, use_cpp = use_cpp)
      cost <- init_ab(cost, spread = spread, min_dist = min_dist, verbose = verbose)
      cost$eps <- eps
      cost
    },
    update = function(cost, Y) {
      D2 <- calc_d2(Y, use_cpp = use_cpp, n_threads = n_threads)
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
cetsne <- function(perplexity, inp_kernel = "gaussian", 
                   symmetrize = "symmetric", eps = .Machine$double.eps,
                   n_threads = 0, use_cpp = FALSE) {
  lreplace(tsne(perplexity, use_cpp = use_cpp, n_threads = n_threads),
           init = function(cost, X, max_iter, verbose = FALSE,
                           ret_extra = c()) {
             cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                              symmetrize = "symmetric", normalize = TRUE,
                              verbose = verbose, ret_extra = ret_extra,
                              n_threads = n_threads, use_cpp = use_cpp)

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


# Bandwidth Experiments --------------------------------------------------

# t-SNE with input kernel bandwidths transferred to output
btsne <- function(perplexity, inp_kernel = "gaussian", beta = NULL,
                  eps = .Machine$double.eps, n_threads = 0, use_cpp = FALSE) {
  lreplace(tsne(perplexity = perplexity, inp_kernel = inp_kernel,
                use_cpp = use_cpp, n_threads = n_threads),
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      ret_extra <- unique(c(ret_extra, 'beta'))
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = "symmetric", normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra,
                       n_threads = n_threads, use_cpp = use_cpp)
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
      W <- calc_d2(Y, use_cpp = use_cpp, n_threads = n_threads)
      W <- 1 / (1 + (cost$beta * W))
      diag(W) <- 0

      cost$Z <- sum(W)
      cost$W <- W
      cost
    }
  )
}

# SSNE with input kernel bandwidths transferred to output
bssne <- function(perplexity, inp_kernel = "gaussian", beta = NULL, 
                  eps = .Machine$double.eps, n_threads = 0, use_cpp = FALSE) {
  lreplace(ssne(perplexity = perplexity, inp_kernel = inp_kernel, 
                use_cpp = use_cpp, n_threads = n_threads),
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
       ret_extra <- unique(c(ret_extra, 'beta'))
       cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                        symmetrize = "symmetric", normalize = TRUE,
                        verbose = verbose, ret_extra = ret_extra,
                        n_threads = n_threads, use_cpp = use_cpp)
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
      cost$Q <- expQ(Y, cost$eps, beta = cost$beta, matrix_normalize = TRUE,
                     use_cpp = use_cpp, n_threads = n_threads)$Q
      cost
    }
  )
}

# ASNE with input kernel bandwidths transferred to output
basne <- function(perplexity, beta = NULL, eps = .Machine$double.eps, 
                  n_threads = 0, use_cpp = FALSE) {
  lreplace(
    asne(perplexity = perplexity, use_cpp = use_cpp, n_threads = n_threads),
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      ret_extra <- unique(c(ret_extra, 'beta'))

      cost <- sne_init(cost, X, perplexity = perplexity,
                       symmetrize = "none", normalize = FALSE,
                       verbose = verbose, ret_extra = ret_extra,
                       n_threads = n_threads, use_cpp = use_cpp)
      # override input bandwidths with fixed beta (although this doesn't do much)
      if (!is.null(beta)) {
        cost$beta <- beta
      }

      cost$eps <- eps
      cost
    },
    update = function(cost, Y) {
      cost$Q <- expQ(Y, eps = cost$eps, beta = cost$beta, is_symmetric = FALSE,
                     use_cpp = use_cpp, n_threads = n_threads)$Q
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
btasne <- function(perplexity, beta = NULL, eps = .Machine$double.eps, 
                   n_threads = 0, use_cpp = FALSE) {
  lreplace(basne(perplexity = perplexity, beta = beta, eps = eps,
                 n_threads = n_threads, use_cpp = use_cpp),
    pfn = kl_cost,
    gr = function(cost, Y) {
    cost <- cost_update(cost, Y)
     W <- cost$W
     cost$G <- k2g(Y, 2 * cost$beta * W * (cost$P - W / cost$Z), 
                   symmetrize = TRUE)
     cost
    },
    update = function(cost, Y) {
      W <- calc_d2(Y, use_cpp = use_cpp, n_threads = n_threads)
      W <- 1 / (1 + cost$beta * W)
      diag(W) <- 0
       
      cost$W <- W
      cost$Z <- rowSums(W)
      cost
    }
  )
}

tasne <- function(perplexity, n_threads = 0, use_cpp = FALSE) {
  lreplace(tsne(perplexity = perplexity, use_cpp = use_cpp,
                n_threads = n_threads),
           init = function(cost, X, max_iter, eps = .Machine$double.eps, 
                           verbose = FALSE, ret_extra = c()) {
             cost <- sne_init(cost, X, perplexity = perplexity,
                              symmetrize = "none", normalize = FALSE,
                              verbose = verbose, ret_extra = ret_extra,
                              n_threads = n_threads, use_cpp = use_cpp)
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
             W <- calc_d2(Y, use_cpp = use_cpp, n_threads = n_threads)
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
tasne <- function(perplexity, inp_kernel = "gaussian", 
                  eps = .Machine$double.eps, n_threads = 0, use_cpp = FALSE) {
  lreplace(tsne(perplexity = perplexity, n_threads = n_threads,
                use_cpp = use_cpp),
  init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
    cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                     symmetrize = "none", normalize = FALSE,
                     verbose = verbose, ret_extra = ret_extra,
                     n_threads = n_threads, use_cpp = use_cpp)
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
    W <- calc_d2(Y, use_cpp = use_cpp, n_threads = n_threads)
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
trmsne <- function(perplexity, inp_kernel = "gaussian", 
                   eps = .Machine$double.eps, n_threads = 0, use_cpp = FALSE) {
  lreplace(
    tsne(perplexity = perplexity, inp_kernel = inp_kernel, 
         n_threads = n_threads, use_cpp = use_cpp),
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = "none", normalize = TRUE, 
                       verbose = verbose, ret_extra = ret_extra, 
                       n_threads = n_threads, use_cpp = use_cpp)
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
tmsne <- function(perplexity, inp_kernel = "gaussian", 
                  eps = .Machine$double.eps, n_threads = 0, use_cpp = FALSE) {
  lreplace(trmsne(perplexity = perplexity, inp_kernel = inp_kernel, 
                  eps = eps, n_threads = n_threads, use_cpp = use_cpp),
           init = function(cost, X, max_iter, verbose = FALSE,
                           ret_extra = c()) {
             cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                              symmetrize = "none", row_normalize = FALSE,
                              normalize = TRUE, verbose = verbose, 
                              ret_extra = ret_extra, n_threads = n_threads,
                              use_cpp = use_cpp)
             cost$eps <- eps
             cost
           }
  )
}

# RSR row-normalize, symmetrize, then row-normalize again
# Might work a tiny bit better than t-ASNE?
trsrsne <- function(perplexity, eps = .Machine$double.eps, n_threads = 0,
                    use_cpp = FALSE) {
  lreplace(tasne(perplexity, use_cpp = use_cpp, n_threads = n_threads),
           init = function(cost, X, max_iter, verbose = FALSE, 
                           ret_extra = c()) {
             cost <- sne_init(cost, X, perplexity = perplexity,
                              symmetrize = "symmetric", normalize = FALSE,
                              verbose = verbose, ret_extra = ret_extra,
                              n_threads = n_threads, use_cpp = use_cpp)
             P <- cost$P
             P <- P / rowSums(P)
             cost$P <- P

             cost$eps <- eps
             cost
           }
  )
}


arsrsne <- function(perplexity, eps = .Machine$double.eps, n_threads = 0,
                    use_cpp = FALSE) {
  lreplace(asne(perplexity, use_cpp = use_cpp, n_threads = n_threads),
           init = function(cost, X, max_iter, verbose = FALSE,
                           ret_extra = c()) {
             cost <- sne_init(cost, X, perplexity = perplexity,
                              symmetrize = "symmetric", normalize = FALSE,
                              verbose = verbose, ret_extra = ret_extra,
                              n_threads = n_threads, use_cpp = use_cpp)
             P <- cost$P
             P <- P / rowSums(P)
             cost$P <- P

             cost$eps <- eps
             cost
           }
  )
}


