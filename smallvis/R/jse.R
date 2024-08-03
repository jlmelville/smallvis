# JSE ---------------------------------------------------------------------

# Lee, J. A., Renard, E., Bernard, G., Dupont, P., & Verleysen, M. (2013).
# Type 1 and 2 mixtures of Kullback-Leibler divergences as cost functions in
# dimensionality reduction based on similarity preservation.
# \emph{Neurocomputing}, \emph{112}, 92-108.
# kappa = 0 behaves like ASNE
# kappa = 1 behaves like NeRV with lambda = 0. Yes that's confusing.
jse <- function(perplexity,
                kappa = 0.5,
                inp_kernel = "gaussian",
                eps = .Machine$double.eps,
                n_threads = 0,
                use_cpp = FALSE) {
  eps0 <- 1e-5
  kappa <- max(kappa, eps0)
  kappa <- min(kappa, 1 - eps0)

  kappa_inv <- 1 / kappa
  m2_kappa_inv <- -2 * kappa_inv
  om_kappa <- 1 - kappa
  om_kappa_inv <- 1 / om_kappa

  lreplace(
    ssne(
      perplexity = perplexity,
      use_cpp = use_cpp,
      n_threads = n_threads
    ),
    init = function(cost,
                    X,
                    max_iter,
                    verbose = FALSE,
                    ret_extra = c()) {
      cost <- sne_init(
        cost,
        X,
        perplexity = perplexity,
        kernel = inp_kernel,
        symmetrize = "none",
        normalize = FALSE,
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
      cost$G <- k2g(Y, K, symmetrize = TRUE)

      cost
    },
    update = function(cost, Y) {
      eps <- cost$eps

      Q <- expQ(
        Y,
        eps = eps,
        is_symmetric = TRUE,
        use_cpp = use_cpp,
        n_threads = n_threads
      )$Q

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


sjse <- function(perplexity,
                 kappa = 0.5,
                 inp_kernel = "gaussian",
                 symmetrize = "symmetric",
                 eps = .Machine$double.eps,
                 n_threads = 0,
                 use_cpp = FALSE) {
  eps0 <- 1e-5
  kappa <- max(kappa, eps0)
  kappa <- min(kappa, 1 - eps0)

  kappa_inv <- 1 / kappa
  m4_kappa_inv <- -4 * kappa_inv
  om_kappa <- 1 - kappa
  om_kappa_inv <- 1 / om_kappa

  lreplace(
    ssne(
      perplexity = perplexity,
      inp_kernel = inp_kernel,
      symmetrize = symmetrize,
      eps = eps,
      n_threads = n_threads,
      use_cpp = use_cpp
    ),
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

      Q <- expQ(
        Y,
        cost$eps,
        is_symmetric = TRUE,
        matrix_normalize = TRUE,
        use_cpp = use_cpp,
        n_threads = n_threads
      )$Q
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

rsrjse <- function(perplexity,
                   kappa = 0.5,
                   eps = .Machine$double.eps,
                   n_threads = 0,
                   use_cpp = FALSE) {
  lreplace(
    jse(
      perplexity = perplexity,
      kappa = kappa,
      use_cpp = use_cpp,
      n_threads = n_threads
    ),
    init = function(cost,
                    X,
                    max_iter,
                    verbose = FALSE,
                    ret_extra = c()) {
      cost <- sne_init(
        cost,
        X,
        perplexity = perplexity,
        symmetrize = "symmetric",
        normalize = FALSE,
        verbose = verbose,
        ret_extra = ret_extra,
        n_threads = n_threads,
        use_cpp = use_cpp
      )
      P <- cost$P
      P <- P / rowSums(P)
      cost$P <- P

      cost$eps <- eps
      cost
    }
  )
}
