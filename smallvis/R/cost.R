# Generic Functions -------------------------------------------------------

logm <- function(m, eps = .Machine$double.eps) {
  diag(m) <- eps
  m <- log(m)
  diag(m) <- 0
  m
}

divm <- function(m, n, eps = .Machine$double.eps) {
  diag(m) <- eps
  m <- m / n
  diag(m) <- 0
  m
}

powm <- function(m, n, eps = .Machine$double.eps) {
  diag(m) <- eps
  m <- m ^ n
  diag(m) <- 0
  m
}

# Convert Force constant to Gradient
k2g <- function(Y, K, symmetrize = FALSE) {
  if (symmetrize) {
    K <- K + t(K)
  }
  Y * colSums(K) - (K %*% Y)
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

largevis <- function(perplexity, inp_kernel = "gaussian", 
                     symmetrize = "symmetric", gamma = 1, gr_eps = 0.1, 
                     normalize = TRUE, eps = 1e-9, row_weight = NULL,
                     use_cpp = FALSE, n_threads = 0) {
  if (!is.null(row_weight)) {
    row_normalize <- row_weight
  }
  else {
    row_normalize <- TRUE
  }
  lreplace(tsne(perplexity, use_cpp = use_cpp, n_threads = n_threads),
     init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
       symmetrize <- match.arg(tolower(symmetrize), true_symmetrize_options())
       cost <- sne_init(cost, X = X, perplexity = perplexity, 
                        symmetrize = symmetrize, kernel = inp_kernel,
                        normalize = normalize, verbose = verbose,
                        row_normalize = row_normalize, ret_extra = ret_extra,
                        n_threads = n_threads, use_cpp = use_cpp)
       cost$eps <- eps
       cost$greps1 <- gr_eps - 1
       cost
     },
     pfn = function(cost, Y) {
       cost <- cost_update(cost, Y)

       P <- cost$P
       eps <- cost$eps
       W <- cost$W

       cost$pcost <- colSums(-P * logm(W, eps) - gamma * log1p(-W + eps))
       cost
     },
     gr = function(cost, Y) {
       cost <- cost_update(cost, Y)
       W <- cost$W
       cost$G <- k2g(Y, 4 * W * (cost$P - ((gamma * W) / (1 + cost$greps1 * W))))
       cost
     },
     update = function(cost, Y) {
       W <- calc_d2(Y, use_cpp = use_cpp, n_threads = n_threads)
       W <- 1 / (1 + W)
       diag(W) <- 0

       cost$W <- W
       cost
     },
     export = cost_export
  )
}


# f-Divergences -----------------------------------------------------------

# Reverse KL divergence
rklsne <- function(perplexity, inp_kernel = "gaussian", 
                   symmetrize = "symmetric", eps = .Machine$double.eps,
                   n_threads = 0, use_cpp = FALSE) {
  lreplace(
    tsne(perplexity = perplexity, use_cpp = use_cpp, n_threads = n_threads),
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      symmetrize <- match.arg(tolower(symmetrize), true_symmetrize_options())
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = symmetrize, normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra,
                       n_threads = n_threads, use_cpp = use_cpp)
      cost$eps <- eps
      cost
    },
    cache_input = function(cost) {
      P <- cost$P
      eps <- cost$eps
      P[P < eps] <- eps
      cost$lP <- logm(P, eps)
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
      W <- calc_d2(Y, use_cpp = use_cpp, n_threads = n_threads)
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
jssne <- function(perplexity, inp_kernel = "gaussian", 
                  symmetrize = "symmetric", eps = .Machine$double.eps,
                  n_threads = 0, use_cpp = FALSE) {
  lreplace(
    tsne(perplexity = perplexity, inp_kernel = inp_kernel, use_cpp == use_cpp,
         n_threads = n_threads),
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      symmetrize <- match.arg(tolower(symmetrize), true_symmetrize_options())
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = symmetrize, normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra,
                       n_threads = n_threads, use_cpp = use_cpp)
      cost$eps <- eps
      cost
    },
    cache_input = function(cost) {
      P <- cost$P
      eps <- cost$eps
      P[P < eps] <- eps
      cost$PlP <- colSums(P * logm(P, eps))
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
      
      W <- calc_d2(Y, use_cpp = use_cpp, n_threads = n_threads)
      W <- 1 / (1 + W)
      diag(W) <- 0
      sumW <- sum(W)
      
      Q <- W / sumW
      Z <- 0.5 * (P + Q)
      
      QlQZ <- logm(Q / Z, eps)
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
chsne <- function(perplexity, inp_kernel = "gaussian", 
                  symmetrize = "symmetric", eps = .Machine$double.eps,
                  n_threads = 0, use_cpp = FALSE) {
  lreplace(
    tsne(perplexity = perplexity, inp_kernel = inp_kernel, use_cpp == use_cpp,
         n_threads = n_threads),
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      symmetrize <- match.arg(tolower(symmetrize), true_symmetrize_options())
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = symmetrize, normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra,
                       n_threads = n_threads, use_cpp = use_cpp)
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
      PQ <- cost$P - cost$Q
      
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
      P <- cost$P
      
      W <- calc_d2(Y, use_cpp = use_cpp, n_threads = n_threads)
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
hlsne <- function(perplexity, inp_kernel = "gaussian", 
                  symmetrize = "symmetric", eps = .Machine$double.eps,
                  n_threads = 0, use_cpp = FALSE) {
  lreplace(
    tsne(perplexity = perplexity, inp_kernel = inp_kernel, use_cpp == use_cpp,
         n_threads = n_threads),
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      symmetrize <- match.arg(tolower(symmetrize), true_symmetrize_options())
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = "symmetric", normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra,
                       n_threads = n_threads, use_cpp = use_cpp)
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
      P <- cost$P
      
      W <- calc_d2(Y, use_cpp = use_cpp, n_threads = n_threads)
      W <- 1 / (1 + W)
      diag(W) <- 0
      Z <- sum(W)
      
      Q <- W / Z
      
      cost$Q <- Q
      cost$Z <- Z
      cost$sQ <- sqrt(Q)
      
      cost
    }
  )
}

# ABSNE -------------------------------------------------------------------

# alpha-beta divergence
absne <- function(perplexity, inp_kernel = "gaussian", 
                  symmetrize = "symmetric", alpha = 1, lambda = 1,
                  eps = .Machine$double.eps, n_threads = 0, use_cpp = FALSE) {
  beta <- lambda - alpha

  eps0 <- 1e-5
  if (abs(alpha) > eps0 && abs(lambda) < eps0) {
    # alpha != 0, beta = -alpha (=> lambda == 0)
    return(absneamb(perplexity = perplexity, inp_kernel = inp_kernel, 
                    symmetrize = symmetrize, alpha = alpha, eps = eps))
  }
  if (abs(alpha) > eps0 && abs(beta) < eps0) {
    # alpha != 0, beta = 0 (=> lambda = alpha)
    return(absneb0(perplexity = perplexity, inp_kernel = inp_kernel, 
                   symmetrize = symmetrize, alpha = alpha, eps = eps))
  }
  if (abs(alpha) < eps0 && abs(beta) > eps0 ) {
    # alpha = 0, beta != 0 (=> lambda = beta)
    return(absnea0(perplexity = perplexity, inp_kernel = inp_kernel, 
                   symmetrize = symmetrize, beta = beta, eps = eps))
  }
  if (abs(alpha) < eps0 && abs(beta) < eps0) {
    # alpha = 0, beta = 0 (=> lambda = 0)
    return(absne00(perplexity = perplexity, inp_kernel = inp_kernel,
                   symmetrize = symmetrize, eps = eps))
  }

  if (abs(lambda) < eps0) {
    lambda <- ifelse(lambda == 0, 1, sign(lambda)) * eps0
  }
  lreplace(
    tsne(perplexity = perplexity, use_cpp = TRUE, n_threads = n_threads),
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      symmetrize <- match.arg(tolower(symmetrize), true_symmetrize_options())
      if (verbose) {
        tsmessage("Using ABSNE with alpha = ", formatC(alpha), 
                  " beta = ", formatC(beta))
      }
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = symmetrize, normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra,
                       n_threads = n_threads, use_cpp = TRUE)
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
      W <- calc_d2(Y, use_cpp = use_cpp, n_threads = n_threads)
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
absneb0 <- function(perplexity, inp_kernel = "gaussian", 
                    symmetrize = "symmetric", alpha = 1,
                    eps = .Machine$double.eps, n_threads = 0, use_cpp = FALSE) {
  lreplace(
    tsne(perplexity = perplexity, use_cpp = use_cpp, n_threads = n_threads),
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      symmetrize <- match.arg(tolower(symmetrize), true_symmetrize_options())
      if (verbose) {
        tsmessage("Using ABSNE with alpha = ", formatC(alpha), 
                  " beta = 0")
      }
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = symmetrize, normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra,
                       n_threads = n_threads, use_cpp = use_cpp)
      
      cost$eps <- eps
      
      cost$inva2 <- 1 / (alpha * alpha)
      cost$invam4 <- 4 / alpha
      cost
    },
    cache_input = function(cost) {
      eps <- cost$eps
      P <- cost$P
      P[P < eps] <- eps
      Pa <- powm(P, alpha, eps)
      cost$Pa <- Pa
      cost$PlPac <- colSums(P * logm(Pa, eps))
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
      W <- calc_d2(Y, use_cpp = use_cpp, n_threads = n_threads)
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
absneamb <- function(perplexity, inp_kernel = "gaussian", 
                     symmetrize = "symmetric", alpha = 1,
                     eps = .Machine$double.eps, n_threads = 0,
                     use_cpp = FALSE) {
  lreplace(
    tsne(perplexity = perplexity, use_cpp = use_cpp, n_threads = n_threads),
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      symmetrize <- match.arg(tolower(symmetrize), true_symmetrize_options())
      if (verbose) {
        tsmessage("Using ABSNE with alpha = ", formatC(alpha), 
                  " beta = -", formatC(alpha))
      }
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = symmetrize, normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra,
                       n_threads = n_threads, use_cpp = use_cpp)

      cost$N <- nrow(cost$P)
      cost$N2 <- (cost$N - 1) * cost$N
      cost$eps <- eps

      cost$inva2 <- 1 / (alpha * alpha)
      cost$invam4 <- 4 / alpha
      cost
    },
    cache_input = function(cost) {
      eps <- cost$eps
      P <- cost$P
      P[P < eps] <- eps
      Pa <- powm(P, alpha, eps)
      cost$Pa <- Pa
      Pa[Pa < eps] <- eps
      cost$lPac <- colSums(logm(Pa, eps))
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
      W <- calc_d2(Y, use_cpp = use_cpp, n_threads = n_threads)
      W <- 1 / (1 + W)
      diag(W) <- 0
      Z <- sum(W)
      Q <- W / Z
      
      cost$Q <- Q
      cost$Z <- Z
      
      Qa <- powm(Q, alpha, cost$eps)
      cost$lQac <- colSums(logm(Qa, cost$eps))
      cost$PadivQa <- divm(cost$Pa, Qa)
      cost$PadivQac <- colSums(cost$PadivQa)
      cost$PadivQas <- sum(cost$PadivQac)
      
      cost
    }
  )
}

# alpha = 0, beta != 0 => lambda = beta
absnea0 <- function(perplexity, inp_kernel = "gaussian", 
                    symmetrize = "symmetric", beta = 1,
                    eps = .Machine$double.eps, n_threads = 0,
                    use_cpp = FALSE) {
  lreplace(
    tsne(perplexity = perplexity, use_cpp = use_cpp, n_threads = n_threads),
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      symmetrize <- match.arg(tolower(symmetrize), true_symmetrize_options())
      if (verbose) {
        tsmessage("Using ABSNE with alpha = 0, beta = ", formatC(beta))
      }
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = symmetrize, normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra,
                       n_threads = n_threads, use_cpp = use_cpp)
      
      cost$eps <- eps
      
      cost$invb2 <- 1 / (beta * beta)
      cost$invbm4 <- 4 / beta
      cost
    },
    cache_input = function(cost) {
      eps <- cost$eps
      P <- cost$P
      P[P < eps] <- eps
      Pb <- powm(P, beta, eps)
      Pb[Pb < eps] <- eps
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
      W <- calc_d2(Y, use_cpp = use_cpp, n_threads = n_threads)
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
absne00 <- function(perplexity, inp_kernel = "gaussian", 
                    symmetrize = "symmetric", eps = .Machine$double.eps,
                    n_threads = 0, use_cpp = FALSE) {
  lreplace(
    tsne(perplexity = perplexity, use_cpp = use_cpp, n_threads = n_threads),
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      symmetrize <- match.arg(tolower(symmetrize), true_symmetrize_options())
      if (verbose) {
        tsmessage("Using ABSNE with alpha = 0, beta = 0")
      }
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = symmetrize, normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra,
                       n_threads = n_threads, use_cpp = use_cpp)
      cost$eps <- eps
      cost
    },
    cache_input = function(cost) {
      eps <- cost$eps
      P <- cost$P
      P[P < eps] <- eps
      lP <- logm(P, eps)
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
      W <- calc_d2(Y, use_cpp = use_cpp, n_threads = n_threads)
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
abssne <- function(perplexity, inp_kernel = "gaussian", 
                   symmetrize = "symmetric", alpha = 1, lambda = 1,
                   eps = .Machine$double.eps, n_threads = 0, use_cpp = FALSE) {
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
    tsne(perplexity = perplexity, use_cpp = use_cpp, n_threads = n_threads),
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      symmetrize <- match.arg(tolower(symmetrize), true_symmetrize_options())
      if (verbose) {
        tsmessage("Using ABSSNE with alpha = ", formatC(alpha), 
                  " beta = ", formatC(beta))
      }
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = symmetrize, normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra,
                       n_threads = n_threads, use_cpp = use_cpp)
      cost$inva4 <- 4 / alpha
      cost$minvab <- -1 / (alpha * beta)
      cost$inval <- 1 / (alpha * lambda)
      cost$ibl <- 1 / (beta * lambda)
      
      cost$eps <- eps
      cost
    },
    cache_input = function(cost) {
      P <- cost$P
      eps <- cost$eps
      P[P < eps] <- eps
      cost$Pa <- powm(P, alpha, eps)
      cost$Plc <- colSums(powm(P, lambda, eps)) * cost$ibl
      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)
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
      eps <- cost$eps
      
      Q <- expQ(Y, cost$eps, is_symmetric = TRUE, matrix_normalize = TRUE,
                use_cpp = use_cpp, n_threads = n_threads)$Q
      
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
gsne <- function(perplexity, lambda = 1, inp_kernel = "gaussian", 
                 symmetrize = "symmetric", eps = .Machine$double.eps,
                 use_cpp = FALSE, n_threads = 0) {
  lreplace(
    tsne(perplexity = perplexity, use_cpp = use_cpp, n_threads = n_threads),
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      symmetrize <- match.arg(tolower(symmetrize), true_symmetrize_options())
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = symmetrize, normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra,
                       n_threads = n_threads, use_cpp = use_cpp)
      cost$eps <- eps
      
      # Phat isn't affected by early exaggeration, so we cache it once only
      if (methods::is(X, "dist")) {
        Phat <- as.matrix(X)
      }
      else {
        Phat <- calc_d2(X, use_cpp = use_cpp, n_threads = n_threads)
      }
      
      Phat <- Phat + 1
      diag(Phat) <- 0
      
      Phat <- Phat / sum(Phat)
      cost$Phat <- Phat
      cost$phlogph <- colSums(Phat * logm(Phat, eps))
      
      cost
    },
    cache_input = function(cost) {
      eps <- cost$eps
      P <- cost$P
      P[P < eps] <- eps
      cost$plogp <- colSums(P * logm(P, eps))

      cost$plamphat <- P - lambda * cost$Phat
      
      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)
      
      eps <- cost$eps
      kl <- cost$plogp - colSums(cost$P * logm(cost$W / cost$Z, eps))
      klhat <- cost$phlogph - colSums(cost$Phat * logm(cost$What / cost$Zhat, eps))
      cost$pcost <- kl + lambda * klhat
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      
      qlampqhat <- (cost$W / cost$Z) - lambda * (cost$What / cost$Zhat)
      cost$G <- k2g(Y, 4 * cost$W * (cost$plamphat - qlampqhat))
      cost
    },
    update = function(cost, Y) {
      W <- calc_d2(Y, use_cpp = use_cpp, n_threads = n_threads)
      What <- 1 + W
      W <- 1 / What
      diag(W) <- 0
      cost$Z <- sum(W)
      cost$W <- W

      diag(What) <- 0
      cost$What <- What
      cost$Zhat <- sum(What)

      cost
    }
  )
}


# Distance Preserving Methods ---------------------------------------------

mmds_init <- function(cost, X, max_iter, eps = .Machine$double.eps, 
                      verbose = FALSE, ret_extra = c(), use_cpp = FALSE,
                      n_threads = 1) {
  tsmessage("Calculating pairwise distances")
  if (methods::is(X, "dist")) {
    cost$R <- as.matrix(X)
  }
  else {
    cost$R <- calc_d(X, use_cpp = use_cpp, n_threads = n_threads)
  }
  cost$eps <- eps
  cost
}

# Metric MDS, minimizing strain.
mmds <- function(eps = .Machine$double.eps, use_cpp = FALSE, n_threads = 1) {
  list(
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      cost <- mmds_init(cost, X, max_iter, eps, verbose, ret_extra, 
                        use_cpp = use_cpp, n_threads = n_threads)
      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)
      cost$pcost <- colSums((cost$R - cost$D) ^ 2)
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      if (use_cpp) {
        cost$G <- mmds_grad_cpp(cost$R, cost$D, Y, eps = eps, n_threads = n_threads)
      }
      else {
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

smmds <- function(eps = .Machine$double.eps, use_cpp = FALSE, n_threads = 1) {
  lreplace(
    mmds(use_cpp = use_cpp, n_threads = n_threads),
    init = function(cost, X, max_iter, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      cost <- mmds_init(cost, X, max_iter, eps, verbose, ret_extra, 
                        use_cpp = use_cpp, n_threads = n_threads)
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
      cost$D2 <- calc_d2(Y, use_cpp = use_cpp, n_threads = n_threads)
      cost
    },
    sentinel = "D2"
  )
}

sammon <- function(eps = .Machine$double.eps, use_cpp = FALSE, n_threads = 1) {
  lreplace(mmds(eps = eps, use_cpp = use_cpp, n_threads = n_threads),
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      cost <- mmds_init(cost, X, max_iter, eps, verbose, ret_extra, 
                        use_cpp = use_cpp, n_threads = n_threads)
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


gmmds <- function(k, eps = .Machine$double.eps, use_cpp = FALSE, n_threads = 0) {
  lreplace(
    mmds(use_cpp = use_cpp, n_threads = n_threads),
    init = function(cost, X, max_iter, verbose = FALSE,
                    ret_extra = c()) {
      cost$R <- geodesic(X, k, n_threads = n_threads, use_cpp = use_cpp,
                         verbose = verbose)
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
ballmmds <- function(f = 0.1, eps = .Machine$double.eps, use_cpp = FALSE,
                     n_threads = 1) {
  lreplace(
    mmds(use_cpp = use_cpp, n_threads = n_threads),
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      cost <- mmds_init(cost = cost, X = X, max_iter = max_iter, eps = eps, 
                        verbose = verbose, ret_extra = ret_extra,
                        use_cpp = use_cpp, n_threads = n_threads)
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
# unless they smaller than the input distance
knnmmds <- function(k, eps = .Machine$double.eps, use_cpp = FALSE, 
                    n_threads = 0) {
  lreplace(
    mmds(use_cpp = use_cpp, n_threads = n_threads),
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      cost <- mmds_init(cost = cost, X = X, max_iter = max_iter, eps = eps, 
                        verbose = verbose, ret_extra = ret_extra,
                        use_cpp = use_cpp, n_threads = n_threads)
      knn <- knn_graph(X = X, k = k, n_threads = n_threads, verbose = verbose)
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

# Elastic Embedding -------------------------------------------------------

# Carreira-Perpinán, M. A. (2010, June).
# The Elastic Embedding Algorithm for Dimensionality Reduction.
# In \emph{Proceedings of the 27th International Conference on Machine Learning (ICML-10)} (pp. 167-174).
# http://faculty.ucmerced.edu/mcarreira-perpinan/papers/icml10.pdf (PDF)
# lambda control the strength of repulsive vs attractive forces
# if neg_weights is true, the repulsive contribution is weighted based on the
# squared input distances. Otherwise, no weighting is applied.
ee <- function(perplexity, lambda = 100, neg_weights = TRUE, 
               inp_kernel = "gaussian", symmetrize = "symmetric", 
               eps = .Machine$double.eps, use_cpp = FALSE, n_threads = 0) {
  list(
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      symmetrize <- match.arg(tolower(symmetrize), true_symmetrize_options())
      if (neg_weights) {
        if (methods::is(X, "dist")) {
          R <- X
        }
        else {
          R <- calc_d(X, use_cpp = use_cpp, n_threads = n_threads)
        }
        cost$Vn <- R / sum(R)
      }
      else {
        cost$Vn <- 1
      }
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = symmetrize, normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra,
                       n_threads = n_threads, use_cpp = use_cpp)
      cost$eps <- eps
      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)

      Vp <- cost$P
      Vn <- cost$Vn
      W <- cost$W
      cost$pcost <- colSums(-Vp * logm(W, cost$eps) + lambda * (Vn * W))
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
      W <- calc_d2(Y, use_cpp = use_cpp, n_threads = n_threads)
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
nerv <- function(perplexity, lambda = 0.9, inp_kernel = "gaussian", 
                 eps = .Machine$double.eps, n_threads = 0, use_cpp = FALSE) {
  lambda2 <- 2 * lambda
  oml <- 1 - lambda
  oml2 <- 2 * oml
  lreplace(
    tsne(perplexity = perplexity, use_cpp = use_cpp, n_threads = n_threads),
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = "none", normalize = FALSE,
                       verbose = verbose, ret_extra = ret_extra,
                       n_threads = n_threads, use_cpp = use_cpp)
      cost$eps <- eps
      cost
    },
    cache_input = function(cost) {
      eps <- cost$eps
      P <- cost$P
      P[P < eps] <- eps
      cost$P <- P
      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)

      P <- cost$P

      kl_fwd <- rowSums(P * cost$lPQ)

      cost$pcost <- lambda * kl_fwd + oml * cost$kl_rev
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      Q <- cost$Q

      # Total K including multiplying by 2 in gradient
      K <- lambda2 * (cost$P - Q) + oml2 * Q * (cost$lPQ + cost$kl_rev)

      cost$G <- k2g(Y, K, symmetrize = TRUE)
      cost
    },
    update = function(cost, Y) {
      eps <- cost$eps

      Q <- expQ(Y, eps, is_symmetric = TRUE, use_cpp = use_cpp, 
                n_threads = n_threads)$Q
      cost$Q <- Q
      
      # Reverse KL gradient
      cost$lPQ <- logm(cost$P / Q, eps)
      # for KLrev we want Q * log(Q/P), so take -ve of log(P/Q)
      cost$kl_rev <- rowSums(Q * -cost$lPQ)
      cost
    },
    sentinel = "Q",
    export = function(cost, val) {
      res <- cost_export(cost, val)
      res
    }
    )
}

snerv <- function(perplexity, lambda = 0.9, inp_kernel = "gaussian", 
                  symmetrize = "symmetric", eps = .Machine$double.eps,
                  n_threads = 0, use_cpp = FALSE) {
  lambda4 <- 4 * lambda
  oml <- 1 - lambda
  oml4 <- 4 * oml
  lreplace(
    ssne(perplexity = perplexity, inp_kernel = inp_kernel, 
         symmetrize = symmetrize, eps = eps, n_threads = n_threads),
    cache_input = function(cost) {
      eps <- cost$eps
      P <- cost$P
      P[P < eps] <- eps
      cost$P <- P
      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)
      
      P <- cost$P
      
      kl_fwd <- colSums(P * cost$lPQ)
      
      cost$pcost <- lambda * kl_fwd + oml * cost$kl_rev
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      Q <- cost$Q
      
      # Total K including multiplying by 4 in gradient
      K <- lambda4 * (cost$P - Q) + oml4 * Q * (cost$lPQ + cost$QlPQs)
      
      cost$G <- k2g(Y, K, symmetrize = FALSE)
      cost
    },
    update = function(cost, Y) {
      eps <- cost$eps
      
      Q <- expQ(Y, eps, is_symmetric = TRUE, matrix_normalize = TRUE,
                use_cpp = use_cpp, n_threads = n_threads)$Q      
      cost$Q <- Q
      
      # Reverse KL gradient
      cost$lPQ <- logm(cost$P / Q, eps)
      cost$kl_rev <- colSums(Q * -cost$lPQ)
      cost$QlPQs <- sum(cost$kl_rev)
      cost
    },
    export = function(cost, val) {
      res <- cost_export(cost, val)
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
jse <- function(perplexity, kappa = 0.5, inp_kernel = "gaussian", 
                eps = .Machine$double.eps, n_threads = 0, use_cpp = FALSE) {
  eps0 <- 1e-5
  kappa <- max(kappa, eps0)
  kappa <- min(kappa, 1 - eps0)

  kappa_inv <- 1 / kappa
  m2_kappa_inv <- -2 * kappa_inv
  om_kappa <- 1 - kappa
  om_kappa_inv <- 1 / om_kappa

  lreplace(
    ssne(perplexity = perplexity, use_cpp = use_cpp, n_threads = n_threads),
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = "none", normalize = FALSE,
                       verbose = verbose, ret_extra = ret_extra,
                       n_threads = n_threads, use_cpp = use_cpp)
      cost$eps <- eps
      cost
    },
    cache_input = function(cost) {
      eps <- cost$eps
      P <- cost$P
      P[P < eps] <- eps
      cost$plogp <- rowSums(P * logm(P, eps))
      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)
      eps <- cost$eps

      cost$pcost <- 
        om_kappa_inv * (cost$plogp - rowSums(cost$P * logm(cost$Z, eps))) +
        kappa_inv * cost$QlQZc
      
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      K <- m2_kappa_inv * cost$Q * (cost$lQZ - cost$QlQZc)
      cost$G <- k2g(Y, K, 
                    symmetrize = TRUE)

      cost
    },
    update = function(cost, Y) {
      eps <- cost$eps

      Q <- expQ(Y, eps = eps, is_symmetric = TRUE, use_cpp = use_cpp,
                n_threads = n_threads)$Q

      Z <- kappa * cost$P + om_kappa * Q
      Z[Z < eps] <- eps
      diag(Z) <- 0
      
      cost$Q <- Q
      cost$Z <- Z
      cost$lQZ <- logm(Q / Z, eps)
      cost$QlQZc <- rowSums(Q * cost$lQZ)

      cost
    },
    export = cost_export
  )
}


sjse <- function(perplexity, kappa = 0.5, inp_kernel = "gaussian", 
                 symmetrize = "symmetric", eps = .Machine$double.eps,
                 n_threads = 0, use_cpp = FALSE) {
  eps0 <- 1e-5
  kappa <- max(kappa, eps0)
  kappa <- min(kappa, 1 - eps0)
  
  kappa_inv <- 1 / kappa
  m4_kappa_inv <- -4 * kappa_inv
  om_kappa <- 1 - kappa
  om_kappa_inv <- 1 / om_kappa
  
  lreplace(
    ssne(perplexity = perplexity, inp_kernel = inp_kernel, 
         symmetrize = symmetrize, eps = eps, n_threads = n_threads,
         use_cpp = use_cpp),
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)
      eps <- cost$eps
      
      cost$pcost <- 
        om_kappa_inv * (cost$plogp - colSums(cost$P * logm(cost$Z, eps))) +
        kappa_inv * cost$QlQZc
      
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      K <- m4_kappa_inv * cost$Q * (cost$lQZ - cost$QlQZs)
      cost$G <- k2g(Y, K, symmetrize = FALSE)

      cost
    },
    update = function(cost, Y) {
      eps <- cost$eps
      
      Q <- expQ(Y, cost$eps, is_symmetric = TRUE, matrix_normalize = TRUE,
                use_cpp = use_cpp, n_threads = n_threads)$Q      
      Z <- kappa * cost$P + om_kappa * Q
      Z[Z < eps] <- eps
      diag(Z) <- 0
      
      cost$Q <- Q
      cost$Z <- Z
      cost$lQZ <- logm(Q / Z, eps)
      cost$QlQZc <- colSums(Q * cost$lQZ)
      cost$QlQZs <- sum(cost$QlQZc)
      
      cost
    },
    export = cost_export
  )
}


rsrnerv <- function(perplexity, lambda = 0.9, eps = .Machine$double.eps,
                    n_threads = 0, use_cpp = FALSE) {
  lreplace(nerv(perplexity = perplexity, lambda = lambda, use_cpp = use_cpp, 
                n_threads = n_threads),
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

rsrjse <- function(perplexity, kappa = 0.5, eps = .Machine$double.eps, 
                   n_threads = 0, use_cpp = FALSE) {
  lreplace(jse(perplexity = perplexity, kappa = kappa, use_cpp = use_cpp,
               n_threads = n_threads),
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

# NeRV with input bandwidths transferred to the output kernel, as in the
# original paper.
bnerv <- function(perplexity, lambda = 0.9, eps = .Machine$double.eps, 
                  n_threads = 0, use_cpp = FALSE) {
  lambda2 <- 2 * lambda
  oml <- 1 - lambda
  oml2 <- 2 * oml
  lreplace(
    nerv(perplexity = perplexity, lambda = lambda, use_cpp = use_cpp, 
         n_threads = n_threads),
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      ret_extra <- unique(c(ret_extra, 'beta'))

      cost <- sne_init(cost, X, perplexity = perplexity,
                       symmetrize = "none", normalize = FALSE,
                       verbose = verbose, ret_extra = ret_extra,
                       n_threads = n_threads, use_cpp = use_cpp)
      cost$eps <- eps
      cost$lambda2b <- lambda2 * cost$beta
      cost$oml2b <- oml2 * cost$beta
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      Q <- cost$Q
      
      # Total K including multiplying by 2 * beta in gradient
      K <- cost$lambda2b * (cost$P - Q) + cost$oml2b * Q * (cost$lPQ + cost$kl_rev)
      
      cost$G <- k2g(Y, K, symmetrize = TRUE)
      cost
      
      cost
    },
    update = function(cost, Y) {
      eps <- cost$eps

      Q <- expQ(Y, eps, beta = cost$beta, is_symmetric = FALSE,
                use_cpp = use_cpp, n_threads = n_threads)$Q
      cost$Q <- Q
      
      # Reverse KL gradient
      cost$lPQ <- logm(cost$P / Q, eps)
      cost$kl_rev <- rowSums(Q * -cost$lPQ)

      cost
    }
  )
}
