# f-Divergences -----------------------------------------------------------

# Reverse KL divergence
rklsne <- function(perplexity,
                   inp_kernel = "gaussian",
                   symmetrize = "symmetric",
                   eps = .Machine$double.eps,
                   n_threads = 0,
                   use_cpp = FALSE) {
  lreplace(
    tsne(
      perplexity = perplexity,
      use_cpp = use_cpp,
      n_threads = n_threads
    ),
    init = function(cost,
                    X,
                    max_iter,
                    verbose = FALSE,
                    ret_extra = c()) {
      symmetrize <- match.arg(tolower(symmetrize), true_symmetrize_options())
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
      cost$G <- k2g(Y, 4 * Q * Q * cost$sumW * (sum(cost$QlQPcs) - cost$lQP))
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
jssne <- function(perplexity,
                  inp_kernel = "gaussian",
                  symmetrize = "symmetric",
                  eps = .Machine$double.eps,
                  n_threads = 0,
                  use_cpp = FALSE) {
  lreplace(
    tsne(
      perplexity = perplexity,
      inp_kernel = inp_kernel,
      use_cpp == use_cpp,
      n_threads = n_threads
    ),
    init = function(cost,
                    X,
                    max_iter,
                    verbose = FALSE,
                    ret_extra = c()) {
      symmetrize <- match.arg(tolower(symmetrize), true_symmetrize_options())
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

      cost$G <- k2g(Y, 2 * Q * Q * cost$sumW * (cost$QlQZs - cost$lQZ))
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
chsne <- function(perplexity,
                  inp_kernel = "gaussian",
                  symmetrize = "symmetric",
                  eps = .Machine$double.eps,
                  n_threads = 0,
                  use_cpp = FALSE) {
  lreplace(
    tsne(
      perplexity = perplexity,
      inp_kernel = inp_kernel,
      use_cpp == use_cpp,
      n_threads = n_threads
    ),
    init = function(cost,
                    X,
                    max_iter,
                    verbose = FALSE,
                    ret_extra = c()) {
      symmetrize <- match.arg(tolower(symmetrize), true_symmetrize_options())
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

      cost$G <- k2g(Y, 4 * Q * Q * Z * (P2Q * invQ - sum(P2Q)))
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
hlsne <- function(perplexity,
                  inp_kernel = "gaussian",
                  symmetrize = "symmetric",
                  eps = .Machine$double.eps,
                  n_threads = 0,
                  use_cpp = FALSE) {
  lreplace(
    tsne(
      perplexity = perplexity,
      inp_kernel = inp_kernel,
      use_cpp == use_cpp,
      n_threads = n_threads
    ),
    init = function(cost,
                    X,
                    max_iter,
                    verbose = FALSE,
                    ret_extra = c()) {
      symmetrize <- match.arg(tolower(symmetrize), true_symmetrize_options())
      cost <- sne_init(
        cost,
        X,
        perplexity = perplexity,
        kernel = inp_kernel,
        symmetrize = "symmetric",
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

      cost$G <- k2g(Y, 4 * Q * Q * Z * (PQ - sPQ))
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
