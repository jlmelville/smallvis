# Optimization algorithms

# Adagrad -----------------------------------------------------------------

adagrad <- function(eta = 0.01, eps = 1e-8, verbose = FALSE) {
  list(
    eta = eta,
    eps = eps,
    init = function(opt, n, k) {
      opt$Ghist <- matrix(0, nrow = n, ncol = k)
      opt
    },
    upd = function(opt, G, iter) {
      Ghist <- opt$Ghist
      eps <- opt$eps
      eta <- opt$eta

      Ghist <- Ghist + G * G
      opt$uY <- -(eta * G) / sqrt(Ghist + eps)

      opt$Ghist <- Ghist
      opt
    }
  )
}

# Adadelta ----------------------------------------------------------------

adadelta <- function(rho = 0.95, eps = 1e-6, verbose = FALSE) {
  list(
    rho = rho,
    eps = eps,
    init = function(opt, n, k) {
      opt$Ghist <- matrix(0, nrow = n, ncol = k)
      opt$uYhist <- matrix(0, nrow = n, ncol = k)

      opt
    },
    upd = function(opt, G, iter) {
      Ghist <- opt$Ghist
      uYhist <- opt$uYhist
      eps <- opt$eps
      rho <- opt$rho

      Ghist <- rho * Ghist + (1 - rho) * G * G
      uY <- -G * sqrt((uYhist + eps) / (Ghist + eps))

      opt$uY <- uY
      opt$Ghist <- Ghist
      opt$uYhist <- rho * uYhist + (1 - rho) * uY * uY
      opt
    }
  )
}

# RMSprop -----------------------------------------------------------------

rmsprop <- function(eta = 0.001, rho = 0.9, eps = 1e-6, verbose = FALSE) {
  list(
    eta = eta,
    rho = rho,
    eps = eps,
    init = function(opt, n, k) {
      opt$Ghist <- matrix(0, nrow = n, ncol = k)

      opt
    },
    upd = function(opt, G, iter) {
      Ghist <- opt$Ghist
      eta <- opt$eta
      eps <- opt$eps
      rho <- opt$rho

      Ghist <- rho * Ghist + (1 - rho) * G * G
      opt$uY <- -G * eta / sqrt(Ghist + eps)

      opt$Ghist <- Ghist
      opt
    }
  )
}

# DBD ---------------------------------------------------------------------

# @param momentum Initial momentum value.
# @param final_momentum Final momentum value.
# @param mom_switch_iter Iteration at which the momentum will switch from
#   \code{momentum} to \code{final_momentum}.
# @param eta Learning rate value.
# @param min_gain Minimum gradient descent step size.
dbd <- function(eta = 500, momentum = 0.5, final_momentum = 0.8,
                mom_switch_iter = 250, min_gain = 0.01, verbose = FALSE) {
  list(
    eta = eta,
    mu = momentum,
    final_momentum = final_momentum,
    mom_switch_iter = mom_switch_iter,
    min_gain = min_gain,
    init = function(opt, n, k) {
      opt$uY <- matrix(0, nrow = n, ncol = k)
      opt$gains <- matrix(1, nrow = n, ncol = k)
      opt
    },
    upd = function(opt, G, iter) {
      gains <- opt$gains
      uY <- opt$uY
      min_gain <- opt$min_gain

      # compare signs of G with -update (== previous G, if we ignore momentum)
      # abs converts TRUE/FALSE to 1/0
      dbd <- abs(sign(G) != sign(uY))
      gains <- (gains + 0.2) * dbd + (gains * 0.8) * (1 - dbd)
      gains[gains < min_gain] <- min_gain
      opt$uY <- opt$mu * opt$uY - opt$eta * gains * G

      opt$gains <- gains

      if (iter == opt$mom_switch_iter && opt$mu != opt$final_momentum) {
        opt$mu <- opt$final_momentum
        if (verbose) {
          message("Iteration #", iter,
                  " switching to final momentum = ", formatC(opt$final_momentum))
        }
      }
      opt
    }
  )
}

# Nesterov DBD ------------------------------------------------------------

ndbd <- function(eta = 500, momentum = 0.5, final_momentum = 0.8,
                 mom_switch_iter = 250, min_gain = 0.01, beta = 1,
                 scale_grad = FALSE,
                 verbose = FALSE) {
  list(
    eta = eta,
    mu = momentum,
    final_momentum = final_momentum,
    mom_switch_iter = mom_switch_iter,
    min_gain = min_gain,
    init = function(opt, n, k) {
      opt$uY <- matrix(0, nrow = n, ncol = k)
      opt$gains <- matrix(1, nrow = n, ncol = k)
      opt
    },
    upd = function(opt, G, iter) {
      gains <- opt$gains
      uY <- opt$uY
      min_gain <- opt$min_gain

      if (scale_grad) {
        G <- G / norm2(G)
      }

      # compare signs of G with -update (== previous G, if we ignore momentum)
      # abs converts TRUE/FALSE to 1/0
      dbd <- abs(sign(G) != sign(uY))
      gains <- (gains + 0.2) * dbd + (gains * 0.8) * (1 - dbd)
      gains[gains < min_gain] <- min_gain

      gd <- opt$eta * gains * G
      old_gd <- opt$old_gd
      # Steepest Descent on first iteration
      if (is.null(old_gd)) {
        uY_new <- -gd
      }
      else {
        # Classical momentum term
        uY_new <- opt$mu * uY - gd
        if (beta != 0) {
          # The momentum bit
          uY_new <- uY_new + beta * opt$mu * (old_gd - gd)
        }
      }

      opt$uY <- uY_new
      opt$gains <- gains
      opt$old_gd <- gd

      if (iter == opt$mom_switch_iter && opt$mu != opt$final_momentum) {
        opt$mu <- opt$final_momentum
        if (verbose) {
          message("Iteration #", iter,
                  " switching to final momentum = ", formatC(opt$final_momentum))
        }
      }
      opt
    }
  )
}

# Adam --------------------------------------------------------------------

adam <- function(eta = 0.002, beta1 = 0.9, beta2 = 0.999, eps = 1e-8,
                 verbose = FALSE) {
  list(
    beta1t = beta1,
    beta2t = beta2,
    eps = eps,
    init = function(opt, n, k) {
      opt$m <- matrix(0, nrow = n, ncol = k)
      opt$v <- matrix(0, nrow = n, ncol = k)

      opt
    },
    upd = function(opt, G, iter) {
      # Uses the slightly simpler formulation from section 2 of the paper
      # which doesn't explicitly store m_hat or v_hat - eps is different
      # between the two versions, though
      m <- opt$m
      v <- opt$v
      beta1t <- opt$beta1t
      beta2t <- opt$beta2t

      # The bias correction part of the update
      bc <- sqrt(1 - beta2t) / (1 - beta1t)

      m <- beta1 * m + (1 - beta1) * G
      v <- beta2 * v + (1 - beta2) * G * G

      gains <- m * bc / (sqrt(v) + opt$eps)

      opt$uY <- -eta * gains
      opt$m <- m
      opt$v <- v
      opt$beta1t <- beta1t * beta1
      opt$beta2t <- beta2t * beta2

      opt
    }
  )
}

# A correction to Adam, with better theoretical convergence properties
# https://openreview.net/forum?id=ryQu7f-RZ
amsgrad <- function(eta = 0.002, beta1 = 0.9, beta2 = 0.999, eps = 1e-8,
                    verbose = FALSE) {
  list(
    beta1t = beta1,
    beta2t = beta2,
    eps = eps,
    init = function(opt, n, k) {
      opt$m <- matrix(0, nrow = n, ncol = k)
      opt$v <- matrix(0, nrow = n, ncol = k)
      opt$v_hat <- matrix(0, nrow = n, ncol = k)

      opt
    },
    upd = function(opt, G, iter) {
      m <- opt$m
      v <- opt$v
      v_hat <- opt$v_hat
      beta1t <- opt$beta1t
      beta2t <- opt$beta2t

      # The bias correction part of the update
      bc <- sqrt(1 - beta2t) / (1 - beta1t)

      m <- beta1 * m + (1 - beta1) * G
      v <- beta2 * v + (1 - beta2) * G * G

      v_hat <- pmax(v_hat, v)

      gains <- m * bc / (sqrt(v_hat) + opt$eps)

      opt$uY <- -eta * gains
      opt$m <- m
      opt$v <- v
      opt$v_hat <- v_hat
      opt$beta1t <- beta1t * beta1
      opt$beta2t <- beta2t * beta2

      opt
    }
  )
}

# AdaMax ------------------------------------------------------------------

adamax <- function(eta = 0.002, beta1 = 0.9, beta2 = 0.999, eps = 1e-8,
                   verbose = FALSE) {
  list(
    beta1 = beta1,
    beta2 = beta2,
    beta1t = beta1,
    eta = eta,
    eps = eps,
    init = function(opt, n, k) {
      opt$m <- matrix(0, nrow = n, ncol = k)
      opt$u <- matrix(0, nrow = n, ncol = k)

      opt
    },
    upd = function(opt, G, iter) {
      m <- opt$m
      u <- opt$u
      eta <- opt$eta
      beta1 <- opt$beta1
      beta2 <- opt$beta2
      beta1t <- opt$beta1t

      m <- beta1 * m + (1 - beta1) * G
      # element-wise max(beta2 * u, |G|)
      u <- apply(cbind(beta2 * u, base::abs(G)), 1, base::max)
      m_hat <- m / (1 - beta1t)

      opt$uY <- -eta * (m_hat / u)
      opt$m <- m
      opt$u <- u
      opt$beta1t <- beta1t * beta1

      opt
    }
  )
}

# Nadam -------------------------------------------------------------------

nadam <- function(eta = 0.002, beta1 = 0.9, beta2 = 0.999, eps = 1e-8,
                  verbose = FALSE) {
  list(
    beta1 = beta1,
    beta2 = beta2,
    beta1t = beta1,
    beta2t = beta2,
    eta = eta,
    eps = eps,
    init = function(opt, n, k) {
      opt$m <- matrix(0, nrow = n, ncol = k)
      opt$v <- matrix(0, nrow = n, ncol = k)

      opt
    },
    upd = function(opt, G, iter) {
      eta <- opt$eta
      eps <- opt$eps
      m <- opt$m
      v <- opt$v
      beta1 <- opt$beta1
      beta2 <- opt$beta2
      beta1t <- opt$beta1t
      beta2t <- opt$beta2t

      m <- beta1 * m + (1 - beta1) * G
      v <- beta2 * v + (1 - beta2) * G * G
      m_hat <- m / (1 - beta1t * beta1)
      v_hat <- v / (1 - beta2t)

      m_hat_nadam <- beta1 * m_hat + (((1 - beta1) * G) / (1 - beta1t))

      opt$uY <- -(eta * m_hat_nadam) / (sqrt(v_hat) + eps)
      opt$m <- m
      opt$v <- v
      opt$beta1t <- beta1t * beta1
      opt$beta2t <- beta2t * beta2

      opt
    }
  )
}

# Steepest Descent --------------------------------------------------------

steepd <- function(eta = 1, verbose = FALSE) {
  list(
    upd = function(opt, G, iter) {
      opt$uY <- -eta * G
      opt
    }
  )
}

# Momentum ------------------------------------------------------

# Classical momentum
mom <- function(eta = 1, mu = 0.9, verbose = FALSE) {
  list(
    upd = function(opt, G, iter) {
      opt$uY <- -eta * G + opt$mu * opt$uY
      opt
    },
    mu = mu,
    uY = 0
  )
}

# Nesterov-style momentum
# if mu is set to a non numeric value, uses the Nesterov schedule
nag <- function(eta = 1, mu = 0.9, verbose = FALSE) {
  list(
    upd = function(opt, G, iter) {
      if (!is.numeric(mu)) {
        mu_t <- 1 - 3 / (5 + iter)
      }
      else {
        mu_t <- mu
      }

      gd <- eta * G
      old_uY <- opt$uY
      old_gd <- opt$old_gd
      if (is.null(opt$old_gd)) {
        uY <- -gd
      }
      else {
        uY <- mu_t * (old_uY + old_gd) - (1 + mu_t) * gd
      }

      opt$old_gd <- gd
      opt$uY <- uY

      opt
    },
    mu = mu
  )
}

# A "generalized" momentum, beta = 0 for classical, beta = 1 for Nesterov
# if mu is set to a non numeric value, uses the Nesterov schedule
gmom <- function(eta = 1, mu = 0.9, beta = 0.5, verbose = FALSE) {
  list(
    upd = function(opt, G, iter) {
      gd <- eta * G
      old_uY <- opt$uY
      old_gd <- opt$old_gd

      if (!is.numeric(mu)) {
        mu_t <- 1 - 3 / (5 + iter)
      }
      else {
        mu_t <- mu
      }

      # Steepest Descent on first iteration
      if (is.null(opt$old_gd)) {
        uY <- -gd
      }
      else {
        # Classical momentum term
        uY <- mu_t * old_uY - gd
        if (beta != 0) {
          uY <- uY + beta * mu_t * (old_gd - gd)
        }
      }

      opt$old_gd <- gd
      opt$uY <- uY

      opt
    },
    mu = mu
  )
}




# Mize Bridge -------------------------------------------------------------

# Adapt the smallvis cost function into the mize form
mizify_cost <- function(cost, nrow) {
  fg <- function(par) {
    dim(par) <- c(nrow, length(par) / nrow)
    cost <- cost_grad(cost, par)
    G <- cost$G

    cost <- cost_point(cost, par)
    pcosts <- cost$pcost
    fn <- sum(pcosts)

    dim(G) <- NULL
    list(fn = fn, gr = G)
  }
  fn <- function(par) {
    fg(par)$fn
  }
  gr <- function(par) {
    fg(par)$gr
  }
  list(
    fn = fn, gr = gr, fg = fg
  )
}

# One step of mize optimization
opt_step_mize <- function(opt, cost_fn, Y, iter) {
  nr <- nrow(Y)
  nc <- ncol(Y)
  fg <- mizify_cost(cost_fn, nr)
  if (nnat(opt$requires_B) && iter == 1) {
    # The Spectral Direction is the positive part of the Hessian, 4L+.
    # mize stores the gradient as a vector, so the Hessian should be a dN x dN
    # block diagonal matrix, with this block repeated d times (d being the
    # output dimension), but shamefully, mize has a hack inside the hessian solve
    # code to recognize when the Hessian is too small for exactly this case.

    # Graph Laplacian , L+ = D+ - W+
    if (is.null(cost_fn$P)) {
      stop("No P matrix in cost for use in approximate Hessian")
    }
    Lp <- diag(colSums(cost_fn$P)) - cost_fn$P
    mu <- min(Lp[Lp > 0]) * 1e-10
    B <- 4 * (Lp + mu)

    fg$hs <- function(par) {
      B
    }
  }
  # Quasi Newton methods need the inverse Hessian
  # so just use the easy-to-invert diagonal
  if (nnat(opt$requires_H) && iter == 1) {
    H <- 4 * colSums(cost_fn$P)
    fg$hi <- function(par) {
      rep(1 / H, nc)
    }
  }
  dim(Y) <- NULL

  if (iter == 1) {
    opt <- mize::mize_init(opt, Y, fg)
  }
  res <- mize::mize_step(opt, Y, fg)
  step_summary <- mize::mize_step_summary(res$opt, res$par, fg, par_old = Y)
  res$opt <- mize::check_mize_convergence(step_summary)

  Y <- res$par
  dim(Y) <- c(nr, nc)

  list(
    Y = Y,
    cost_fn = cost_fn,
    opt = res$opt,
    f = res$f,
    G = res$g
  )
}

# Generic Functions -------------------------------------------------------

# Creates the optimizer: mize or internal
opt_create <- function(optlist, verbose = FALSE) {
  name <- optlist[[1]]
  optlist[[1]] <- NULL

  if (tolower(name) %in% c("adagrad", "adadelta", "rmsprop", "dbd",
                           "ndbd", "adam", "adamax", "nadam", "amsgrad",
                           "steepd", "mom", "nag", "gmom")) {
    optlist$verbose <- verbose
    opt <- do.call(get(name), optlist)
    opt$smallvis_step <- opt_step_internal
  }
  else {
    if (tolower(name) == "specd") {
      optlist$method <- "PHESS"
    }
    else {
      optlist$method <- name
    }
    optlist$max_iter <- Inf
    approx_hess <- optlist$approx_hess
    optlist$approx_hess <- NULL
    opt <- do.call(mize::make_mize, optlist)
    opt$name <- name
    opt$smallvis_step <- opt_step_mize
    if (tolower(name) %in% c("specd")) {
      opt$requires_B <- TRUE
    }
    # Use the diagonal of the Spectral Direction Hessian approximation
    # to initialize the quasi-Newton methods
    if (nnat(approx_hess) && tolower(name) %in% c("bfgs", "l-bfgs", "sr1")) {
      opt$requires_H <- TRUE
    }
  }

  opt
}

# Initializes the optimizer - called on the first iteration
opt_init <- function(opt, n, ndim) {
  if (!is.null(opt$init)) {
    opt <- opt$init(opt, n, ndim)
  }
  opt
}

# Runs one "step" of optimization - could be multiple f/g evaluations with
# mize optimizers
opt_step <- function(opt, cost_fn, Y, iter) {
  opt$smallvis_step(opt, cost_fn, Y, iter)
}

# One step of optimization using the non-mize optimizers
opt_step_internal <- function(opt, cost_fn, Y, iter) {
  if (iter == 1) {
    opt <- opt_init(opt, nrow(Y), ncol(Y))
  }

  cost_fn <- cost_grad(cost_fn, Y)

  opt <- opt$upd(opt, cost_fn$G, iter)
  cost_fn <- cost_clear(cost_fn)
  list(
    Y = Y + opt$uY,
    cost_fn = cost_fn,
    opt = opt,
    G = cost_fn$G
  )
}

