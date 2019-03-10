# Generic Functions -------------------------------------------------------

logm <- function(m, eps = .Machine$double.eps) {
  diag(m) <- eps
  m <- log(m)
  diag(m) <- 0
  m
}

divm <- function(m, n, eps = .Machine$double.eps) {
  m <- m / n
  diag(m) <- 0
  m
}

powm <- function(m, n, eps = .Machine$double.eps) {
  m <- m ^ n
  diag(m) <- 0
  m
}

# Convert Force constant to Gradient
k2g <- function(Y, K, symmetrize = FALSE) {
  if (symmetrize) {
    K <- K + t(K)
  }
  Y * rowSums(K) - (K %*% Y)
}

cost_init <- function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
  if (!is.null(cost$init)) {
    cost <- cost$init(cost, X, verbose = verbose, ret_extra = ret_extra,
                      max_iter = max_iter)
  }
  cost_cache_input(cost)
}

cost_grad <- function(cost, Y) {
  cost$gr(cost, Y)
}

cost_point <- function(cost, Y) {
  cost$pfn(cost, Y)
}

# Clear values dependent on Y (or other parameters)
# This is the mechanism by which Gradient and Function evaluations will
# detect that they need recalculate distances, weights, probabilities etc.
cost_clear <- function(cost) {
  if (!is.null(cost$sentinel)) {
    cost[[cost$sentinel]] <- NULL
  }
  else if (is.null(cost$clear)) {
    cost <- cost$clear(cost)
  }
  cost
}

# Update matrices dependent on Y
# Check for a NULL value of the sentinel to avoid unnecessary recalculation
cost_update <- function(cost, Y) {
  if (!is.null(cost$sentinel)) {
    if (is.null(cost[[cost$sentinel]])) {
      cost <- cost$update(cost, Y)
    }
  }
  else {
    cost <- cost$update(cost, Y)
  }
  cost
}

# Calculate the cost function value. If opt_res$f is non-NULL, use that
cost_eval <- function(cost, Y, opt_res = NULL) {
  if (is.null(opt_res) || is.null(opt_res$f)) {
    cost <- cost_point(cost, Y)
    pcosts <- cost$pcost
    cost_val <- sum(pcosts)
  }
  else {
    cost_val <- opt_res$f
  }

  list(
    cost = cost,
    value = cost_val
  )
}

# Default export of values associated with a method
# If the value is associated with the cost list, e.g. cost$P, then asking
# for val = "P" or val = "p" will return it. Otherwise, returns NULL
cost_export <- function(cost, val) {
  res <- NULL
  if (!is.null(cost[[val]])) {
    res <- cost[[val]]
  }
  else if (!is.null(cost[[toupper(val)]])) {
    res <- cost[[toupper(val)]]
  }
  res
}

cost_cache_input <- function(cost) {
  if (!is.null(cost$cache_input)) {
    cost <- cost$cache_input(cost)
  }
  cost  
}

start_exaggerating <- function(cost, exaggeration_factor) {
  if (!is.null(cost$exaggerate) && exaggeration_factor != 1) {
    cost <- cost$exaggerate(cost, exaggeration_factor)
  }
  cost_cache_input(cost)
}

stop_exaggerating <- function(cost, exaggeration_factor) {
  if (!is.null(cost$exaggerate) && exaggeration_factor != 1) {
    cost <- cost$exaggerate(cost, 1 / exaggeration_factor)
  }
  cost_cache_input(cost)
}

# LargeVis ----------------------------------------------------------------

# NB This version doesn't normalize the input P, despite what the paper
# indicates (source code of the current implementation doesn't seem to either)
largevis <- function(perplexity, gamma = 7, gr_eps = 0.1) {
  lreplace(tsne(perplexity),
     init = function(cost, X, max_iter, eps = 1e-9, verbose = FALSE, ret_extra = c()) {
       cost <- sne_init(cost, X = X, perplexity = perplexity, symmetrize = "symmetric",
                        normalize = TRUE, verbose = verbose,
                        ret_extra = ret_extra)
       cost$eps <- eps
       cost
     },
     pfn = function(cost, Y) {
       cost <- cost_update(cost, Y)

       P <- cost$P
       eps <- cost$eps
       W <- cost$W

       cost$pcost <- colSums(-P * log(W + eps) - gamma * log1p(-W + eps))
       cost
     },
     gr = function(cost, Y) {
       cost <- cost_update(cost, Y)

       W <- cost$W
       cost$G <- k2g(Y, 4 * W * (cost$P - ((gamma * W) / ((1 - W) + gr_eps))))
       cost
     },
     update = function(cost, Y) {
       W <- dist2(Y)
       W <- 1 / (1 + W)
       diag(W) <- 0

       cost$W <- W
       cost
     },
     export = cost_export
  )
}

# UMAP --------------------------------------------------------------------

umap <- function(perplexity, spread = 1, min_dist = 0.001, gr_eps = 0.1) {
  lreplace(tsne(perplexity),
    init = function(cost, X, max_iter, eps = 1e-9, verbose = FALSE, ret_extra = c()) {
      ab_params <- find_ab_params(spread = spread, min_dist = min_dist)
      a <- ab_params[1]
      b <- ab_params[2]
      if (verbose) {
        message("Umap curve parameters = ", formatC(a), ", ", formatC(b))
      }
      cost$a <- a
      cost$b <- b

      P <- smooth_knn_distances(X, k = perplexity, tol = 1e-5,
                                verbose = verbose)$P
      P <- fuzzy_set_union(P)
      cost$P <- P
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
      D2 <- dist2(Y)
      D2[D2 < 0] <- 0

      W <- 1 / (1 + cost$a * D2 ^ cost$b)
      diag(W) <- 0

      cost$W <- W
      cost$D2 <- D2
      cost
    },
    export = cost_export
  )
}

# UMAP with the output kernel fixed to the t-distribution
tumap <- function(perplexity, gr_eps = 0.1) {
  lreplace(umap(perplexity),
    init = function(cost, X, max_iter, eps = 1e-9, verbose = FALSE, ret_extra = c()) {
      cost$eps <- eps

      P <- smooth_knn_distances(X, k = perplexity, tol = 1e-5,
                                verbose = verbose)$P
      P <- fuzzy_set_union(P)
      cost$P <- P
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)

      cost$G <- k2g(Y, 4 * (cost$W / ((1 - cost$W) + cost$eps + gr_eps)) * (cost$P - cost$W))
      cost
    },
    update = function(cost, Y) {
      W <- dist2(Y)
      W <- 1 / (1 + W)
      diag(W) <- 0

      cost$W <- W
      cost
    }
  )
}

# t-UMAP where output and input affinities are normalized
ntumap <- function(perplexity, gr_eps = 0.1) {
  lreplace(tumap(perplexity),
    init = function(cost, X, max_iter, eps = 1e-9, verbose = FALSE, ret_extra = c()) {
      cost$eps <- eps

      P <- smooth_knn_distances(X, k = perplexity, tol = 1e-5,
                                 verbose = verbose)$P
      P <- fuzzy_set_union(P)

      P <- P / sum(P)
      cost$P <- P

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

# f-Divergences -----------------------------------------------------------

# Reverse KL divergence
rklsne <- function(perplexity, inp_kernel = "gaussian") {
  lreplace(
    tsne(perplexity = perplexity, inp_kernel = inp_kernel),
    init = function(cost, X, max_iter, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = "symmetric", normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra)
      cost$eps <- eps
      cost
    },
    cache_input = function(cost) {
      cost$lP <- logm(cost$P, cost$eps)
      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)
      cost$pcost <- cost$QlQPcs
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      Q <- cost$Q
      cost$G <- k2g(Y,  4 * Q * Q * cost$sumW * (sum(cost$QlQPcs) - cost$lQP))
      cost
    },
    update = function(cost, Y) {
      W <- dist2(Y)
      W <- 1 / (1 + W)
      diag(W) <- 0
      sumW <- sum(W)
      Q <- W / sumW

      cost$lQP <- logm(Q, cost$eps) - cost$lP
      cost$QlQPcs <- colSums(Q * cost$lQP)
      cost$Q <- Q
      cost$sumW <- sumW
      
      cost
    }
  )
}

# Jensen-Shannon divergence
jssne <- function(perplexity, inp_kernel = "gaussian") {
  lreplace(
    tsne(perplexity = perplexity, inp_kernel = inp_kernel),
    init = function(cost, X, max_iter, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = "symmetric", normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra)
      cost$eps <- eps
      cost
    },
    cache_input = function(cost) {
      P <- cost$P
      cost$PlP <- colSums(P * logm(P, cost$eps))
      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)
      cost$pcost <- 0.5 * (cost$PlP + cost$QlQZ - colSums(cost$P * logm(cost$Z, cost$eps)))
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      Q <- cost$Q

      cost$G <- k2g(Y,  2 * Q * Q * cost$sumW * (cost$QlQZs - cost$lQZ))
      cost
    },
    update = function(cost, Y) {
      eps <- cost$eps
      P <- cost$P
      
      W <- dist2(Y)
      W <- 1 / (1 + W)
      diag(W) <- 0
      sumW <- sum(W)
      
      Q <- W / sumW
      Z <- 0.5 * (P + Q)
      
      QlQZ <- logm(Q / Z)
      cost$QlQZ <- colSums(Q * QlQZ)
      cost$QlQZs <- sum(cost$QlQZ)
      
      cost$lQZ <- QlQZ
      cost$sumW <- sumW
      cost$Q <- Q
      cost$Z <- Z 
      
      cost
    }
  )
}

# Chi-squared divergence
chsne <- function(perplexity, inp_kernel = "gaussian") {
  lreplace(
    tsne(perplexity = perplexity, inp_kernel = inp_kernel),
    init = function(cost, X, max_iter, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = "symmetric", normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra)
      cost$eps <- eps
      cost
    },
    cache_input = function(cost) {
      P <- cost$P
      cost$P2 <- P * P
      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)
      P <- cost$P
      Q <- cost$Q      

      PQ <- P - Q
      
      cost$pcost <- colSums(PQ * PQ * cost$invQ)
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      
      Q <- cost$Q
      Z <- cost$Z
      invQ <- cost$invQ
      
      P2Q <- cost$P2 * invQ
      
      cost$G <- k2g(Y,  4 * Q * Q * Z * (P2Q * invQ - sum(P2Q)))
      cost
    },
    update = function(cost, Y) {
      eps <- cost$eps
      P <- cost$P
      
      W <- dist2(Y)
      W <- 1 / (1 + W)
      diag(W) <- 0
      Z <- sum(W)
      
      Q <- W / Z

      cost$Q <- Q
      cost$Z <- Z
      invQ <- 1 / Q
      diag(invQ) <- 0
      cost$invQ <- invQ
      
      cost
    }
  )
}

# Hellinger distance divergence
hlsne <- function(perplexity, inp_kernel = "gaussian") {
  lreplace(
    tsne(perplexity = perplexity, inp_kernel = inp_kernel),
    init = function(cost, X, max_iter, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = "symmetric", normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra)
      cost$eps <- eps
      cost
    },
    cache_input = function(cost) {
      cost$sP <- sqrt(cost$P)
      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)

      PQ <- cost$sP - cost$sQ
      
      cost$pcost <- colSums(PQ * PQ)
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      sP <- cost$sP
      Q <- cost$Q
      Z <- cost$Z
      sQ <- cost$sQ
      
      sPQ <- sum(sP * sQ)
      PQ <- sP / sQ
      diag(PQ) <- 0
      
      cost$G <- k2g(Y,  4 * Q * Q * Z * (PQ - sPQ))
      cost
    },
    update = function(cost, Y) {
      eps <- cost$eps
      P <- cost$P
      
      W <- dist2(Y)
      W <- 1 / (1 + W)
      diag(W) <- 0
      Z <- sum(W)
      
      Q <- W / Z
      
      cost$Q <- Q
      cost$Z <- Z
      cost$sQ <- sqrt(Q)
      invQ <- 1 / Q
      diag(invQ) <- 0
      cost$invQ <- invQ
      
      cost
    }
  )
}

# ABSNE -------------------------------------------------------------------

# alpha-beta divergence
absne <- function(perplexity, inp_kernel = "gaussian", alpha = 1, lambda = 1) {
  beta <- lambda - alpha

  eps0 <- 1e-5
  if (abs(alpha) > eps0 && abs(lambda) < eps0) {
    # alpha != 0, beta = -alpha (=> lambda == 0)
    return(absneamb(perplexity = perplexity,
                    inp_kernel = inp_kernel, alpha = alpha))
  }
  if (abs(alpha) > eps0 && abs(beta) < eps0) {
    # alpha != 0, beta = 0 (=> lambda = alpha)
    return(absneb0(perplexity = perplexity,
                    inp_kernel = inp_kernel, alpha = alpha))
  }
  if (abs(alpha) < eps0 && abs(beta) > eps0 ) {
    # alpha = 0, beta != 0 (=> lambda = beta)
    return(absnea0(perplexity = perplexity,
                   inp_kernel = inp_kernel, beta = beta))
  }
  if (abs(alpha) < eps0 && abs(beta) < eps0) {
    # alpha = 0, beta = 0 (=> lambda = 0)
    return(absne00(perplexity = perplexity, inp_kernel = inp_kernel))
  }

  if (abs(lambda) < eps0) {
    lambda <- ifelse(lambda == 0, 1, sign(lambda)) * eps0
  }
  lreplace(
    tsne(perplexity = perplexity, inp_kernel = inp_kernel),
    init = function(cost, X, max_iter, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      if (verbose) {
        tsmessage("Using ABSNE with alpha = ", formatC(alpha), 
                  " beta = ", formatC(beta))
      }
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = "symmetric", normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra)
      cost$inva4 <- 4 / alpha
      cost$minvab <- -1 / (alpha * beta)
      cost$inval <- 1 / (alpha * lambda)

      cost$eps <- eps
      cost
    },
    cache_input = function(cost) {
      P <- cost$P
      eps <- cost$eps
      cost$Pa <- powm(P, alpha, eps)
      cost$Plc <- colSums(powm(P, lambda, eps)) / (beta * lambda)
      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)
      cost$pcost <- cost$minvab * cost$PaQbc + cost$Plc + cost$inval * cost$Qlc
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      Q <- cost$Q
      cost$G <- k2g(Y, cost$inva4 * cost$Z * Q *  
                      (cost$PaQb - cost$Ql + Q * (cost$Qls - cost$PaQbs)))
      cost
    },
    update = function(cost, Y) {
      W <- dist2(Y)
      W <- 1 / (1 + W)
      diag(W) <- 0
      Z <- sum(W)
      Q <- W / Z

      eps <- cost$eps
      cost$PaQb <- cost$Pa * powm(Q, beta, eps)
      cost$PaQbc <- colSums(cost$PaQb)
      cost$PaQbs <- sum(cost$PaQbc)

      cost$Q <- Q
      cost$Ql <- powm(Q, lambda, eps)
      cost$Qlc <- colSums(cost$Ql)
      cost$Qls <- sum(cost$Qlc)  
          
      cost$Z <- Z
      cost
    }
  )
}

# alpha != 0, beta = 0 => lambda = alpha
absneb0 <- function(perplexity, inp_kernel = "gaussian", alpha = 1) {
  lreplace(
    tsne(perplexity = perplexity, inp_kernel = inp_kernel),
    init = function(cost, X, max_iter, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      if (verbose) {
        tsmessage("Using ABSNE with alpha = ", formatC(alpha), 
                  " beta = 0")
      }
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = "symmetric", normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra)
      
      cost$eps <- eps
      
      cost$inva2 <- 1 / (alpha * alpha)
      cost$invam4 <- 4 / alpha
      cost
    },
    cache_input = function(cost) {
      eps <- cost$eps
      Pa <- powm(cost$P, alpha, eps)
      cost$Pa <- Pa
      cost$PlPac <- colSums(cost$P * logm(Pa, eps))
      cost$Pac <- colSums(Pa)
      cost$Pas <- sum(cost$Pac)
      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)
      cost$pcost <- cost$inva2 * 
        (cost$PlPac - colSums(cost$Pa * logm(cost$Qa, cost$eps)) - 
           cost$Pac + cost$Qac)
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      Q <- cost$Q
      
      cost$G <- k2g(Y, cost$invam4 * Q * cost$Z *
                      (cost$Pa - cost$Qa + Q * (cost$Qas - cost$Pas)))
      cost
    },
    update = function(cost, Y) {
      W <- dist2(Y)
      W <- 1 / (1 + W)
      diag(W) <- 0
      Z <- sum(W)
      Q <- W / Z
      
      cost$Q <- Q
      cost$Z <- Z
      
      Qa <- powm(Q, alpha, cost$eps)
      cost$Q <- Q
      cost$Qa <- Qa
      cost$Qac <- colSums(Qa)
      cost$Qas <- sum(cost$Qac)
      
      cost
    }
  )
}

# alpha = -beta != 0 => lambda = 0
absneamb <- function(perplexity, inp_kernel = "gaussian", alpha = 1) {
  lreplace(
    tsne(perplexity = perplexity, inp_kernel = inp_kernel),
    init = function(cost, X, max_iter, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      if (verbose) {
        tsmessage("Using ABSNE with alpha = ", formatC(alpha), 
                  " beta = -", formatC(alpha))
      }
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = "symmetric", normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra)

      cost$N <- nrow(cost$P)
      cost$N2 <- (cost$N - 1) * cost$N
      cost$eps <- eps

      cost$inva2 <- 1 / (alpha * alpha)
      cost$invam4 <- 4 / alpha
      cost
    },
    cache_input = function(cost) {
      eps <- cost$eps
      cost$Pa <- powm(cost$P, alpha, eps)
      cost$lPac <- colSums(logm(cost$Pa, eps))
      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)
      cost$pcost <- cost$inva2 * (cost$lQac - cost$lPac + cost$PadivQac - cost$N)
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      Q <- cost$Q

      cost$G <- k2g(Y, cost$invam4 * Q * cost$Z *
                    (cost$PadivQa - 1 + Q * (cost$N2 - cost$PadivQas)))
      cost
    },
    update = function(cost, Y) {
      W <- dist2(Y)
      W <- 1 / (1 + W)
      diag(W) <- 0
      Z <- sum(W)
      Q <- W / Z
      
      cost$Q <- Q
      cost$Z <- Z
      
      Qa <- Q ^ alpha
      cost$lQac <- colSums(logm(Qa, cost$eps))
      cost$PadivQa <- divm(cost$Pa, Qa)
      cost$PadivQac <- colSums(cost$PadivQa)
      cost$PadivQas <- sum(cost$PadivQac)
      
      cost
    }
  )
}

# alpha = 0, beta != 0 => lambda = beta
absnea0 <- function(perplexity, inp_kernel = "gaussian", beta = 1) {
  lreplace(
    tsne(perplexity = perplexity, inp_kernel = inp_kernel),
    init = function(cost, X, max_iter, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      if (verbose) {
        tsmessage("Using ABSNE with alpha = 0, beta = ", formatC(beta))
      }
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = "symmetric", normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra)
      
      cost$eps <- eps
      
      cost$invb2 <- 1 / (beta * beta)
      cost$invbm4 <- 4 / beta
      cost
    },
    cache_input = function(cost) {
      eps <- cost$eps
      Pb <- powm(cost$P, beta, eps)
      cost$Pbc <- colSums(Pb)
      cost$lPb <- logm(Pb, eps)
      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)
      Qb <- cost$Qb
      cost$pcost <- cost$invb2 * 
        (cost$QblQbc - cost$QblPbc + cost$Pbc - cost$Qbc)
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      Q <- cost$Q
      
      cost$G <- k2g(Y, cost$invbm4 * Q * cost$Z *
                      (cost$QblPb - cost$QblQb + Q * (cost$QblQbs - cost$QblPbs)))
      cost
    },
    update = function(cost, Y) {
      W <- dist2(Y)
      W <- 1 / (1 + W)
      diag(W) <- 0
      Z <- sum(W)
      Q <- W / Z
      
      cost$Q <- Q
      cost$Z <- Z
      
      Qb <- powm(Q, beta, cost$eps)
      cost$Qbc <- colSums(Qb)
      
      cost$QblPb <- Qb * cost$lPb
      cost$QblPbc <- colSums(cost$QblPb)
      cost$QblPbs <- sum(cost$QblPb)
      
      cost$QblQb <- Qb * logm(Qb, cost$eps)
      cost$QblQbc <- colSums(cost$QblQb)
      cost$QblQbs <- sum(cost$QblQb)
      
      cost
    }
  )
}

# alpha = 0, beta = 0 => lambda = 0
absne00 <- function(perplexity, inp_kernel = "gaussian") {
  lreplace(
    tsne(perplexity = perplexity, inp_kernel = inp_kernel),
    init = function(cost, X, max_iter, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      if (verbose) {
        tsmessage("Using ABSNE with alpha = 0, beta = 0")
      }
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = "symmetric", normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra)
      
      cost$eps <- eps
      cost
    },
    cache_input = function(cost) {
      eps <- cost$eps
      lP <- logm(cost$P, eps)
      cost$lP <- lP
      cost$lPs <- sum(lP)
      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)
      pcost <- cost$lP - cost$lQ
      cost$pcost <- 0.5 * colSums(pcost * pcost)
      
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      Q <- cost$Q
      
      cost$G <- k2g(Y, 4 * Q * cost$Z *
                      (cost$lP - cost$lQ + Q * (cost$lQs - cost$lPs)))
      cost
    },
    update = function(cost, Y) {
      W <- dist2(Y)
      W <- 1 / (1 + W)
      diag(W) <- 0
      Z <- sum(W)
      Q <- W / Z
      
      cost$Q <- Q
      cost$Z <- Z
      
      lQ <- logm(Q, cost$eps)
      cost$lQ <- lQ
      cost$lQs <- sum(lQ)
      
      cost
    }
  )
}

# alpha-beta divergence
abssne <- function(perplexity, inp_kernel = "gaussian", alpha = 1, lambda = 1) {
  beta <- lambda - alpha
  
  eps0 <- 1e-5

  if (abs(alpha) < eps0) {
    alpha <- ifelse(alpha == 0, 1, sign(alpha)) * eps0
  }
  if (abs(beta) < eps0) {
    beta <- ifelse(beta == 0, 1, sign(beta)) * eps0
  }
  if (abs(lambda) < eps0) {
    lambda <- ifelse(lambda == 0, 1, sign(lambda)) * eps0
  }
  lreplace(
    tsne(perplexity = perplexity, inp_kernel = inp_kernel),
    init = function(cost, X, max_iter, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      if (verbose) {
        tsmessage("Using ABSSNE with alpha = ", formatC(alpha), 
                  " beta = ", formatC(beta))
      }
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = "symmetric", normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra)
      cost$inva4 <- 4 / alpha
      cost$minvab <- -1 / (alpha * beta)
      cost$inval <- 1 / (alpha * lambda)
      
      cost$eps <- eps
      cost
    },
    cache_input = function(cost) {
      P <- cost$P
      eps <- cost$eps
      cost$Pa <- powm(P, alpha, eps)
      cost$Plc <- colSums(powm(P, lambda, eps)) / (beta * lambda)
      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)
      # browser()
      cost$pcost <- cost$minvab * cost$PaQbc + cost$Plc + cost$inval * cost$Qlc
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      cost$G <- k2g(Y, cost$inva4 *  
                      (cost$PaQb - cost$Ql + cost$Q * (cost$Qls - cost$PaQbs)))
      cost
    },
    update = function(cost, Y) {
      W <- dist2(Y)
      W <- exp(-W)
      diag(W) <- 0
      Q <- W / sum(W)
      
      eps <- cost$eps
      cost$PaQb <- cost$Pa * powm(Q, beta, eps)
      cost$PaQbc <- colSums(cost$PaQb)
      cost$PaQbs <- sum(cost$PaQbc)
      
      cost$Q <- Q
      cost$Ql <- powm(Q, lambda, eps)
      cost$Qlc <- colSums(cost$Ql)
      cost$Qls <- sum(cost$Qlc)  
      
      cost
    }
  )
}

# Other Divergences -------------------------------------------------------

# global-SNE
gsne <- function(perplexity, lambda = 1, inp_kernel = "gaussian") {
  lreplace(
    tsne(perplexity = perplexity, inp_kernel = inp_kernel),
    init = function(cost, X, max_iter, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = "symmetric", normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra)
      cost$eps <- eps
      
      # Phat isn't affected by early exaggeration, so we cache it once only
      if (methods::is(X, "dist")) {
        Phat <- as.matrix(X)
      }
      else {
        Phat <- safe_dist2(X)
      }
      
      Phat <- Phat + 1
      diag(Phat) <- 0
      
      Phat <- Phat / sum(Phat)
      cost$Phat <- Phat
      cost$phlogph <- colSums(Phat * log(Phat + eps))
      
      cost
    },
    cache_input = function(cost) {
      eps <- cost$eps
      P <- cost$P
      cost$plogp <- colSums(P * log((P + eps)))

      cost$plamphat <- P - lambda * cost$Phat
      
      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)
      
      eps <- cost$eps
      
      P <- cost$P
      invZ <- cost$invZ
      W <- cost$W
      kl <- cost$plogp - colSums(P * log((W * invZ) + eps))
      
      Phat <- cost$Phat
      invZhat <- cost$invZhat
      What <- cost$What
      klhat <- cost$phlogph - colSums(Phat * log((What * invZhat) + eps))
      
      cost$pcost <- kl + lambda * klhat
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      
      qlampqhat <- (cost$W * cost$invZ) - lambda * (cost$What * cost$invZhat)
      cost$G <- k2g(Y, 4 * cost$W * (cost$plamphat - qlampqhat))
      cost
    },
    update = function(cost, Y) {
      W <- dist2(Y)
      What <- 1 + W
      W <- 1 / What
      diag(W) <- 0
      cost$invZ <- 1 / sum(W)
      cost$W <- W

      diag(What) <- 0
      cost$What <- What
      cost$invZhat <- 1 / sum(What)

      cost
    }
  )
}


# Distance Preserving Methods ---------------------------------------------

mmds_init <- function(cost, X, max_iter, eps = .Machine$double.eps, verbose = FALSE,
                      ret_extra = c()) {
  if (methods::is(X, "dist")) {
    cost$R <- as.matrix(X)
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
      cost <- cost_update(cost, Y)
      cost$pcost <- colSums((cost$R - cost$D) ^ 2)
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      cost$G <- k2g(Y, -4 * (cost$R - cost$D) / (cost$D + cost$eps))
      cost
    },
    update = function(cost, Y) {
      cost$D <- sqrt(safe_dist2(Y))
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

smmds <- function() {
  lreplace(
    mmds(),
    init = function(cost, X, max_iter, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      cost <- mmds_init(cost, X, max_iter, eps, verbose, ret_extra)
      cost$R2 <- cost$R * cost$R
      cost$R <- NULL
      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)
      cost$pcost <- colSums((cost$R2 - cost$D2) ^ 2)
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
      cost$D2 <- dist2(Y)
      cost
    },
    sentinel = "D2"
  )
}

sammon <- function() {
  lreplace(mmds(),
    init = function(cost, X, max_iter, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      cost <- mmds_init(cost, X, max_iter, eps, verbose, ret_extra)
      cost$rsum_inv <- 1 / sum(cost$R)
      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)
      cost$pcost <- colSums((cost$R - cost$D) ^ 2 / (cost$R + cost$eps)) * cost$rsum_inv
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      cost$G <- k2g(Y, -4 * cost$rsum_inv * (cost$R - cost$D) / (cost$R * cost$D + cost$eps))
      cost
    }
  )
}


gmmds <- function(k) {
  lreplace(
    mmds(),
    init = function(cost, X, max_iter, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      cost$R <- geodesic(X, k, verbose = verbose)

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
ballmmds <- function(f = 0.1) {
  lreplace(
    mmds(),
    init = function(cost, X, max_iter, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      cost <- mmds_init(cost = cost, X = X, max_iter = max_iter, eps = eps, verbose = verbose,
                        ret_extra = ret_extra)

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
      cost$D <- sqrt(safe_dist2(Y))
      cost
    }
  )
}

# Create the symmetrized knn graph, don't correct non-neighborhood distances
# unless they smaller than the input distance
knnmmds <- function(k) {
  lreplace(
    mmds(),
    init = function(cost, X, max_iter, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      cost <- mmds_init(cost = cost, X = X, max_iter = max_iter, eps = eps, verbose = verbose,
                        ret_extra = ret_extra)
      knn <- knn_graph(X = X, k = k)
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
      cost$D <- sqrt(safe_dist2(Y))
      cost
    }
  )
}

# Elastic Embedding -------------------------------------------------------

# Carreira-Perpinán, M. A. (2010, June).
# The Elastic Embedding Algorithm for Dimensionality Reduction.
# In \emph{Proceedings of the 27th International Conference on Machine Learning (ICML-10)} (pp. 167-174).
# http://faculty.ucmerced.edu/mcarreira-perpinan/papers/icml10.pdf (PDF)
# lambda control the strength of repulsive vs attractive forces
# if neg_weights is true, the repulsive contribution is weighted based on the
# squared input distances. Otherwise, no weighting is applied.
ee <- function(perplexity, lambda = 100, neg_weights = TRUE) {
  list(
    init = function(cost, X, max_iter, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {

      if (neg_weights) {
        if (methods::is(X, "dist")) {
          R <- X
        }
        else {
          R <- sqrt(safe_dist2(X))
        }
        cost$Vn <- R / sum(R)
      }
      else {
        cost$Vn <- 1
      }
      cost <- sne_init(cost, X, perplexity = perplexity,
                       symmetrize = "symmetric", normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra)
      cost$eps <- eps
      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)

      Vp <- cost$P
      Vn <- cost$Vn
      W <- cost$W
      eps <- cost$eps
      cost$pcost <- colSums(-Vp * log(W + eps) + lambda * (Vn * W))
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)

      Vp <- cost$P
      Vn <- cost$Vn
      cost$G <- k2g(Y,  4 * (Vp - lambda * Vn * cost$W))
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
    },
    update = function(cost, Y) {
      W <- dist2(Y)
      W <- exp(-W)
      diag(W) <- 0
      cost$W <- W
      cost
    },
    sentinel = "W"
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
    init = function(cost, X, max_iter, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      cost <- sne_init(cost, X, perplexity = perplexity,
                       symmetrize = "none", normalize = FALSE,
                       verbose = verbose, ret_extra = ret_extra)
      cost$eps <- eps
      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)

      P <- cost$P

      kl_fwd <- rowSums(P * cost$lPQ)

      cost$pcost <- lambda * kl_fwd + (1 - lambda) * cost$kl_rev
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)

      eps <- cost$eps
      P <- cost$P
      Q <- cost$Q

      # Forward KL gradient
      K <- lambda * (P - Q)
      # Total K
      K <- K + (1 - lambda) * (Q * (cost$lPQ + cost$kl_rev))

      cost$G <- k2g(Y, 2 * K, symmetrize = TRUE)
      cost
    },
    update = function(cost, Y) {
      eps <- cost$eps

      W <- dist2(Y)
      W <- exp(-W)
      # Particularly for low lambda, need to make sure W is never all zero
      W[W < eps] <- eps
      diag(W) <- 0
      invZ <- 1 / colSums(W)
      Q <- W * invZ

      # Reverse KL gradient
      lPQ <- log((cost$P + eps) / (Q + eps))
      # for KLrev we want Q * log(Q/P), so take -ve of log(P/Q)
      kl_rev <- rowSums(Q * -lPQ)

      cost$invZ <- invZ
      cost$Q <- Q
      cost$kl_rev <- kl_rev
      cost$lPQ <- lPQ

      cost
    },
    sentinel = "Q",
    export = function(cost, val) {
      res <- cost_export(cost, val)

      if (is.null(res)) {
        switch(val,
               w = {
                 res <- cost$Q / cost$invZ
               })
      }
      res
    }
    )
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
    init = function(cost, X, max_iter, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      cost <- sne_init(cost, X, perplexity = perplexity,
                       symmetrize = "none", normalize = FALSE,
                       verbose = verbose, ret_extra = ret_extra)
      cost$eps <- eps
      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)

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
      cost <- cost_update(cost, Y)

      K <- kappa_inv * (cost$Q * (cost$lZQ + cost$kl_rev))
      cost$G <- k2g(Y, 2 * K, symmetrize = TRUE)
      cost
    },
    update = function(cost, Y) {
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

      cost$Q <- Q
      cost$Z <- Z
      cost$lZQ <- lZQ
      cost$kl_rev <- kl_rev

      cost
    },
    export = cost_export
  )
}


rsrnerv <- function(perplexity, lambda = 0.9) {
  lreplace(nerv(perplexity = perplexity, lambda = lambda),
           init = function(cost, X, max_iter, eps = .Machine$double.eps, verbose = FALSE,
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

rsrjse <- function(perplexity, kappa = 0.5) {
  lreplace(jse(perplexity = perplexity, kappa = kappa),
           init = function(cost, X, max_iter, eps = .Machine$double.eps, verbose = FALSE,
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

# NeRV with input bandwidths transferred to the output kernel, as in the
# original paper.
bnerv <- function(perplexity, lambda = 0.9) {
  lreplace(
    nerv(perplexity = perplexity, lambda = lambda),
    init = function(cost, X, max_iter, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      ret_extra <- unique(c(ret_extra, 'beta'))

      cost <- sne_init(cost, X, perplexity = perplexity,
                       symmetrize = "none", normalize = FALSE,
                       verbose = verbose, ret_extra = ret_extra)
      cost$eps <- eps
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)

      eps <- cost$eps
      Q <- cost$Q

      # Forward KL gradient
      K <- lambda * (cost$P - Q)
      # Total K
      K <- K + (1 - lambda) * (Q * (cost$lPQ + cost$kl_rev))

      cost$G <- k2g(Y, 2 * cost$beta * K, symmetrize = TRUE)
      cost
    },
    update = function(cost, Y) {
      eps <- cost$eps

      W <- dist2(Y)
      W <- exp(-W * cost$beta)

      W[W < eps] <- eps
      diag(W) <- 0
      invZ <- 1 / rowSums(W)
      Q <- W * invZ

      # Reverse KL gradient
      lPQ <- log((cost$P + eps) / (Q + eps))
      kl_rev <- rowSums(Q * -lPQ)

      cost$invZ <- invZ
      cost$Q <- Q
      cost$kl_rev <- kl_rev
      cost$lPQ <- lPQ

      cost
    })
}
