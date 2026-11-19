# ABSNE -------------------------------------------------------------------

absne_limit <- function(alpha, beta, lambda) {
  threshold <- 1e-5
  if (abs(alpha) <= threshold) {
    if (abs(beta) <= threshold) {
      return("joint")
    }
    return("alpha")
  }
  if (abs(lambda) <= threshold) {
    return("lambda")
  }
  if (abs(beta) <= threshold) {
    return("beta")
  }
  "regular"
}

# alpha-beta divergence
absne <- function(
  perplexity,
  inp_kernel = "gaussian",
  symmetrize = "symmetric",
  alpha = 1,
  lambda = 1,
  eps = .Machine$double.eps,
  n_threads = 0,
  use_cpp = FALSE
) {
  beta <- lambda - alpha
  limit <- absne_limit(alpha, beta, lambda)

  if (limit == "lambda") {
    # alpha != 0, beta = -alpha (=> lambda == 0)
    return(
      absneamb(
        perplexity = perplexity,
        inp_kernel = inp_kernel,
        symmetrize = symmetrize,
        alpha = alpha,
        eps = eps,
        n_threads = n_threads,
        use_cpp = use_cpp
      )
    )
  }
  if (limit == "beta") {
    # alpha != 0, beta = 0 (=> lambda = alpha)
    return(
      absneb0(
        perplexity = perplexity,
        inp_kernel = inp_kernel,
        symmetrize = symmetrize,
        alpha = alpha,
        eps = eps,
        n_threads = n_threads,
        use_cpp = use_cpp
      )
    )
  }
  if (limit == "alpha") {
    # alpha = 0, beta != 0 (=> lambda = beta)
    return(
      absnea0(
        perplexity = perplexity,
        inp_kernel = inp_kernel,
        symmetrize = symmetrize,
        beta = beta,
        eps = eps,
        n_threads = n_threads,
        use_cpp = use_cpp
      )
    )
  }
  if (limit == "joint") {
    # alpha = 0, beta = 0 (=> lambda = 0)
    return(
      absne00(
        perplexity = perplexity,
        inp_kernel = inp_kernel,
        symmetrize = symmetrize,
        eps = eps,
        n_threads = n_threads,
        use_cpp = use_cpp
      )
    )
  }

  lreplace(
    tsne(
      perplexity = perplexity,
      use_cpp = use_cpp,
      n_threads = n_threads
    ),
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      symmetrize <- match.arg(tolower(symmetrize), true_symmetrize_options())
      if (verbose) {
        tsmessage(
          "Using ABSNE with alpha = ",
          formatC(alpha),
          " beta = ",
          formatC(beta)
        )
      }
      cost <- sne_init(
        cost,
        X,
        perplexity = perplexity,
        kernel = inp_kernel,
        symmetrize = symmetrize,
        normalize = TRUE,
        verbose = verbose,
        ret_extra = ret_extra,
        n_threads = n_threads,
        use_cpp = use_cpp
      )
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
      cost$G <- k2g(
        Y,
        cost$inva4 *
          cost$Z *
          Q *
          (cost$PaQb - cost$Ql + Q * (cost$Qls - cost$PaQbs))
      )
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
absneb0 <- function(
  perplexity,
  inp_kernel = "gaussian",
  symmetrize = "symmetric",
  alpha = 1,
  eps = .Machine$double.eps,
  n_threads = 0,
  use_cpp = FALSE
) {
  lreplace(
    tsne(
      perplexity = perplexity,
      use_cpp = use_cpp,
      n_threads = n_threads
    ),
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      symmetrize <- match.arg(tolower(symmetrize), true_symmetrize_options())
      if (verbose) {
        tsmessage("Using ABSNE with alpha = ", formatC(alpha), " beta = 0")
      }
      cost <- sne_init(
        cost,
        X,
        perplexity = perplexity,
        kernel = inp_kernel,
        symmetrize = symmetrize,
        normalize = TRUE,
        verbose = verbose,
        ret_extra = ret_extra,
        n_threads = n_threads,
        use_cpp = use_cpp
      )

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
      cost$PlPac <- colSums(Pa * logm(Pa, eps))
      cost$Pac <- colSums(Pa)
      cost$Pas <- sum(cost$Pac)
      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)
      cost$pcost <- cost$inva2 *
        (cost$PlPac -
          colSums(cost$Pa * logm(cost$Qa, cost$eps)) -
          cost$Pac +
          cost$Qac)
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      Q <- cost$Q

      cost$G <- k2g(
        Y,
        cost$invam4 *
          Q *
          cost$Z *
          (cost$Pa - cost$Qa + Q * (cost$Qas - cost$Pas))
      )
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
absneamb <- function(
  perplexity,
  inp_kernel = "gaussian",
  symmetrize = "symmetric",
  alpha = 1,
  eps = .Machine$double.eps,
  n_threads = 0,
  use_cpp = FALSE
) {
  lreplace(
    tsne(
      perplexity = perplexity,
      use_cpp = use_cpp,
      n_threads = n_threads
    ),
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      symmetrize <- match.arg(tolower(symmetrize), true_symmetrize_options())
      if (verbose) {
        tsmessage(
          "Using ABSNE with alpha = ",
          formatC(alpha),
          " beta = -",
          formatC(alpha)
        )
      }
      cost <- sne_init(
        cost,
        X,
        perplexity = perplexity,
        kernel = inp_kernel,
        symmetrize = symmetrize,
        normalize = TRUE,
        verbose = verbose,
        ret_extra = ret_extra,
        n_threads = n_threads,
        use_cpp = use_cpp
      )

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
      cost$pcost <- cost$inva2 *
        (cost$lQac - cost$lPac + cost$PadivQac - (cost$N - 1))
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      Q <- cost$Q

      cost$G <- k2g(
        Y,
        cost$invam4 *
          Q *
          cost$Z *
          (cost$PadivQa - 1 + Q * (cost$N2 - cost$PadivQas))
      )
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
absnea0 <- function(
  perplexity,
  inp_kernel = "gaussian",
  symmetrize = "symmetric",
  beta = 1,
  eps = .Machine$double.eps,
  n_threads = 0,
  use_cpp = FALSE
) {
  lreplace(
    tsne(
      perplexity = perplexity,
      use_cpp = use_cpp,
      n_threads = n_threads
    ),
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      symmetrize <- match.arg(tolower(symmetrize), true_symmetrize_options())
      if (verbose) {
        tsmessage("Using ABSNE with alpha = 0, beta = ", formatC(beta))
      }
      cost <- sne_init(
        cost,
        X,
        perplexity = perplexity,
        kernel = inp_kernel,
        symmetrize = symmetrize,
        normalize = TRUE,
        verbose = verbose,
        ret_extra = ret_extra,
        n_threads = n_threads,
        use_cpp = use_cpp
      )

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
      cost$pcost <- cost$invb2 *
        (cost$QblQbc - cost$QblPbc + cost$Pbc - cost$Qbc)
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      Q <- cost$Q

      cost$G <- k2g(
        Y,
        cost$invbm4 *
          Q *
          cost$Z *
          (cost$QblPb - cost$QblQb + Q * (cost$QblQbs - cost$QblPbs))
      )
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
absne00 <- function(
  perplexity,
  inp_kernel = "gaussian",
  symmetrize = "symmetric",
  eps = .Machine$double.eps,
  n_threads = 0,
  use_cpp = FALSE
) {
  lreplace(
    tsne(
      perplexity = perplexity,
      use_cpp = use_cpp,
      n_threads = n_threads
    ),
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      symmetrize <- match.arg(tolower(symmetrize), true_symmetrize_options())
      if (verbose) {
        tsmessage("Using ABSNE with alpha = 0, beta = 0")
      }
      cost <- sne_init(
        cost,
        X,
        perplexity = perplexity,
        kernel = inp_kernel,
        symmetrize = symmetrize,
        normalize = TRUE,
        verbose = verbose,
        ret_extra = ret_extra,
        n_threads = n_threads,
        use_cpp = use_cpp
      )
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

      cost$G <- k2g(
        Y,
        4 * Q * cost$Z * (cost$lP - cost$lQ + Q * (cost$lQs - cost$lPs))
      )
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
abssne <- function(
  perplexity,
  inp_kernel = "gaussian",
  symmetrize = "symmetric",
  alpha = 1,
  lambda = 1,
  eps = .Machine$double.eps,
  n_threads = 0,
  use_cpp = FALSE
) {
  beta <- lambda - alpha
  limit <- absne_limit(alpha, beta, lambda)
  exponent <- switch(
    limit,
    joint = 0,
    alpha = beta,
    beta = alpha,
    lambda = 0,
    regular = lambda
  )
  lreplace(
    tsne(
      perplexity = perplexity,
      use_cpp = use_cpp,
      n_threads = n_threads
    ),
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      symmetrize <- match.arg(tolower(symmetrize), true_symmetrize_options())
      if (verbose) {
        tsmessage(
          "Using ABSSNE with alpha = ",
          formatC(alpha),
          " beta = ",
          formatC(beta)
        )
      }
      cost <- sne_init(
        cost,
        X,
        perplexity = perplexity,
        kernel = inp_kernel,
        symmetrize = symmetrize,
        normalize = TRUE,
        verbose = verbose,
        ret_extra = ret_extra,
        n_threads = n_threads,
        use_cpp = use_cpp
      )
      cost$eps <- eps
      cost
    },
    cache_input = function(cost) {
      P <- pmax(cost$P, cost$eps)
      if (limit == "regular") {
        cost$Pa <- powm(P, alpha, cost$eps)
        cost$P_power_sum <- colSums(powm(P, exponent, cost$eps))
      } else {
        cost$logP <- logm(P, cost$eps)
        if (exponent != 0) {
          cost$P_power <- powm(P, exponent, cost$eps)
        }
      }
      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)
      if (limit == "regular") {
        cost$pcost <- (alpha /
          lambda *
          cost$P_power_sum +
          beta / lambda * colSums(cost$Q_power) -
          colSums(cost$cross_power)) /
          (alpha * beta)
        return(cost)
      }
      ratio <- cost$log_ratio
      divergence <- switch(
        limit,
        joint = 0.5 * ratio^2,
        alpha = (cost$P_power -
          cost$Q_power -
          beta * weighted_log_ratio(cost$Q_power, ratio)) /
          beta^2,
        beta = (cost$Q_power -
          cost$P_power +
          alpha * weighted_log_ratio(cost$P_power, ratio)) /
          alpha^2,
        lambda = {
          r <- alpha * ratio
          value <- expm1(r) - r
          value[is.infinite(r)] <- Inf
          value / alpha^2
        }
      )
      cost$pcost <- colSums(divergence)
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      ratio <- cost$log_ratio
      force <- switch(
        limit,
        joint = ratio,
        alpha = weighted_log_ratio(cost$Q_power, ratio),
        beta = (cost$P_power - cost$Q_power) / alpha,
        lambda = expm1(alpha * ratio) / alpha,
        regular = (cost$cross_power - cost$Q_power) / alpha
      )
      # The second term differentiates the Gaussian probabilities' normalizer.
      cost$G <- k2g(Y, 4 * (force - cost$Q * sum(force)))
      cost
    },
    update = function(cost, Y) {
      Q <- expQ(
        Y,
        cost$eps,
        is_symmetric = TRUE,
        matrix_normalize = TRUE,
        use_cpp = use_cpp,
        n_threads = n_threads
      )$Q

      if (exponent != 0) {
        cost$Q_power <- powm(Q, exponent, cost$eps)
      }
      if (limit == "regular") {
        cost$cross_power <- cost$Pa * powm(Q, beta, cost$eps)
      } else {
        cost$log_ratio <- cost$logP - logm(Q, cost$eps)
      }
      cost$Q <- Q
      cost
    },
    sentinel = "Q"
  )
}

weighted_log_ratio <- function(weights, ratio) {
  result <- weights * ratio
  # A zero weight contributes zero even if its log ratio is infinite.
  result[weights == 0] <- 0
  result
}
