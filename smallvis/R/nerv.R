# NeRV --------------------------------------------------------------------

# Venna, J., Peltonen, J., Nybo, K., Aidos, H., & Kaski, S. (2010).
# Information retrieval perspective to nonlinear dimensionality reduction for
# data visualization.
# \emph{Journal of Machine Learning Research}, \emph{11}, 451-490.
#
# Unlike original publication, won't transfer input precisions to output kernel
# lambda = 1 gives ASNE results
# default lambda = 0.9 from "Majorization-Minimization for Manifold Embedding"
# Yang, Peltonen, Kaski 2015
nerv <- function(perplexity,
                 lambda = 0.9,
                 inp_kernel = "gaussian",
                 eps = .Machine$double.eps,
                 n_threads = 0,
                 use_cpp = FALSE) {
  lambda2 <- 2 * lambda
  oml <- 1 - lambda
  oml2 <- 2 * oml
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

      Q <- expQ(
        Y,
        eps,
        is_symmetric = TRUE,
        use_cpp = use_cpp,
        n_threads = n_threads
      )$Q
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

snerv <- function(perplexity,
                  lambda = 0.9,
                  inp_kernel = "gaussian",
                  symmetrize = "symmetric",
                  eps = .Machine$double.eps,
                  n_threads = 0,
                  use_cpp = FALSE) {
  lambda4 <- 4 * lambda
  oml <- 1 - lambda
  oml4 <- 4 * oml
  lreplace(
    ssne(
      perplexity = perplexity,
      inp_kernel = inp_kernel,
      symmetrize = symmetrize,
      eps = eps,
      n_threads = n_threads
    ),
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

      Q <- expQ(
        Y,
        eps,
        is_symmetric = TRUE,
        matrix_normalize = TRUE,
        use_cpp = use_cpp,
        n_threads = n_threads
      )$Q
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

# NeRV with input bandwidths transferred to the output kernel, as in the
# original paper.
bnerv <- function(perplexity,
                  lambda = 0.9,
                  eps = .Machine$double.eps,
                  n_threads = 0,
                  use_cpp = FALSE) {
  lambda2 <- 2 * lambda
  oml <- 1 - lambda
  oml2 <- 2 * oml
  lreplace(
    nerv(
      perplexity = perplexity,
      lambda = lambda,
      use_cpp = use_cpp,
      n_threads = n_threads
    ),
    init = function(cost,
                    X,
                    max_iter,
                    verbose = FALSE,
                    ret_extra = c()) {
      ret_extra <- unique(c(ret_extra, "beta"))

      cost <- sne_init(
        cost,
        X,
        perplexity = perplexity,
        symmetrize = "none",
        normalize = FALSE,
        verbose = verbose,
        ret_extra = ret_extra,
        n_threads = n_threads,
        use_cpp = use_cpp
      )
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

      Q <- expQ(
        Y,
        eps,
        beta = cost$beta,
        is_symmetric = FALSE,
        use_cpp = use_cpp,
        n_threads = n_threads
      )$Q
      cost$Q <- Q

      # Reverse KL gradient
      cost$lPQ <- logm(cost$P / Q, eps)
      cost$kl_rev <- rowSums(Q * -cost$lPQ)

      cost
    }
  )
}

rsrnerv <- function(perplexity,
                    lambda = 0.9,
                    eps = .Machine$double.eps,
                    n_threads = 0,
                    use_cpp = FALSE) {
  lreplace(
    nerv(
      perplexity = perplexity,
      lambda = lambda,
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
