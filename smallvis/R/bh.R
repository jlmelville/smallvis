bhhssne <- function(
  perplexity,
  alpha = 0.5,
  inp_kernel = "nngaussian",
  symmetrize = "symmetric",
  eps = .Machine$double.eps,
  n_threads = 0,
  theta = 0.5,
  nn = "approximate"
) {
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha)) {
    stop("Barnes-Hut HSSNE requires a finite scalar alpha")
  }
  bh_sne(
    perplexity,
    inp_kernel,
    symmetrize,
    eps,
    n_threads,
    theta,
    nn,
    alpha = max(alpha, 1e-8),
    normalize = TRUE,
    lambda = 0
  )
}

bhtee <- function(
  perplexity,
  inp_kernel = "nngaussian",
  symmetrize = "symmetric",
  lambda = 0.01,
  eps = .Machine$double.eps,
  n_threads = 0,
  theta = 0.5,
  nn = "approximate"
) {
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda)) {
    stop("Barnes-Hut TEE requires a finite scalar lambda")
  }
  bh_sne(
    perplexity,
    inp_kernel,
    symmetrize,
    eps,
    n_threads,
    theta,
    nn,
    alpha = 1,
    normalize = FALSE,
    lambda = lambda
  )
}

bh_sne <- function(
  perplexity,
  inp_kernel,
  symmetrize,
  eps,
  n_threads,
  theta,
  nn,
  alpha,
  normalize,
  lambda
) {
  inp_kernel <- match.arg(tolower(inp_kernel), c("nngaussian", "knn"))
  symmetrize <- match.arg(tolower(symmetrize), true_symmetrize_options())
  nn <- match.arg(tolower(nn), c("exact", "approximate"))
  if (
    !is.numeric(n_threads) ||
      length(n_threads) != 1L ||
      !is.finite(n_threads) ||
      n_threads < 0 ||
      n_threads > .Machine$integer.max ||
      n_threads != trunc(n_threads)
  ) {
    stop("n_threads must be a nonnegative integer")
  }
  list(
    ndim = 2,
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      cost <- sne_init(
        cost,
        X,
        perplexity = perplexity,
        kernel = inp_kernel,
        symmetrize = symmetrize,
        normalize = normalize,
        verbose = verbose,
        ret_extra = ret_extra,
        n_threads = n_threads,
        nn = nn,
        sparse = TRUE
      )
      P <- Matrix::drop0(methods::as(cost$P, "generalMatrix"))
      cost$P <- P
      cost$invN <- if (normalize) 1 else 1 / sum(P)
      cost
    },
    cache_input = function(cost) {
      # TEE keeps the input-only term at the original affinity scale.
      if (!normalize && !is.null(cost$constant)) {
        return(cost)
      }
      # The dense objectives floor the input-only terms. Count implicit
      # zeros analytically so that flooring never fills the sparse matrix.
      P <- cost$P
      n <- nrow(P)
      values <- pmax(P@x, eps)
      missing <- n * (n - 1.0) - length(values)
      constant <- sum(values * log(values))
      if (eps > 0) {
        constant <- constant + missing * eps * log(eps)
      }
      if (!normalize) {
        constant <- constant - lambda * (sum(values) + (missing + n) * eps)
      }
      cost$constant <- constant * cost$invN
      cost
    },
    pfn = function(cost, Y) {
      P <- cost$P
      score <- bh_sne_score_cpp(
        P@i,
        P@p,
        P@x,
        Y,
        theta,
        0,
        n_threads,
        alpha,
        normalize,
        lambda
      )
      cost$pcost <- cost$constant - cost$invN * score
      cost
    },
    gr = function(cost, Y) {
      P <- cost$P
      G <- bh_sne_gradient_cpp(
        P@i,
        P@p,
        P@x,
        Y,
        theta,
        0,
        n_threads,
        alpha,
        normalize,
        lambda
      )
      cost$G <- (4 * cost$invN) * G
      cost
    },
    export = cost_export,
    update = function(cost, Y) cost,
    exaggerate = function(cost, exaggeration_factor) {
      cost$P <- cost$P * exaggeration_factor
      cost
    }
  )
}
