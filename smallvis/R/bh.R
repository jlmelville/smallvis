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
  validate_tee_lambda(lambda)
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

bhlargevis <- function(
  perplexity,
  inp_kernel = "nngaussian",
  symmetrize = "symmetric",
  gamma = 1,
  gr_eps = 0.1,
  normalize = TRUE,
  eps = 1e-9,
  row_weight = NULL,
  n_threads = 0,
  theta = 0.5,
  nn = "approximate"
) {
  if (!is.numeric(gamma) || length(gamma) != 1L || !is.finite(gamma)) {
    stop("Barnes-Hut LargeVis requires a finite scalar gamma")
  }
  if (
    !is.numeric(gr_eps) ||
      length(gr_eps) != 1L ||
      !is.finite(gr_eps) ||
      gr_eps < 0
  ) {
    stop("Barnes-Hut LargeVis requires a finite nonnegative gr_eps")
  }
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps < 0) {
    stop("Barnes-Hut LargeVis requires a finite nonnegative eps")
  }
  list(
    ndim = 2,
    gradient_only = gr_eps != 0,
    init = bh_sparse_init(
      perplexity,
      inp_kernel,
      symmetrize,
      normalize,
      n_threads,
      nn,
      row_normalize = if (is.null(row_weight)) TRUE else row_weight
    ),
    pfn = function(cost, Y) {
      P <- cost$P
      cost$pcost <- -bh_largevis_score_cpp(
        P@i,
        P@p,
        P@x,
        Y,
        theta,
        n_threads,
        gamma,
        eps
      )
      cost
    },
    gr = function(cost, Y) {
      P <- cost$P
      cost$G <- 4 *
        bh_largevis_gradient_cpp(
          P@i,
          P@p,
          P@x,
          Y,
          theta,
          n_threads,
          gamma,
          gr_eps
        )
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

bhtumap <- function(
  perplexity,
  inp_kernel = "nngaussian",
  symmetrize = "umap",
  gr_eps = 0.1,
  eps = 1e-9,
  row_weight = NULL,
  n_threads = 0,
  theta = 0.5,
  nn = "approximate"
) {
  if (
    !is.numeric(gr_eps) ||
      length(gr_eps) != 1L ||
      !is.finite(gr_eps) ||
      gr_eps < 0
  ) {
    stop("Barnes-Hut t-UMAP requires a finite nonnegative gr_eps")
  }
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("Barnes-Hut t-UMAP requires a finite positive eps")
  }
  list(
    ndim = 2,
    gradient_only = gr_eps != 0,
    init = bh_sparse_init(
      perplexity,
      inp_kernel,
      symmetrize,
      normalize = FALSE,
      n_threads = n_threads,
      nn = nn,
      row_normalize = if (is.null(row_weight)) FALSE else row_weight,
      kernels = c("nngaussian", "knn", "skd")
    ),
    cache_input = function(cost) {
      p <- cost$P@x
      # At P = W = 0 the two complement terms cancel.
      cost$Cp <- sum(
        p * log(p + eps) + (1 - p) * log1p(eps - p) - log1p(eps)
      )
      cost
    },
    pfn = function(cost, Y) {
      P <- cost$P
      cost$pcost <- cost$Cp -
        bh_tumap_score_cpp(
          P@i,
          P@p,
          P@x,
          Y,
          theta,
          n_threads,
          eps
        )
      cost
    },
    gr = function(cost, Y) {
      P <- cost$P
      cost$G <- 4 *
        bh_tumap_gradient_cpp(
          P@i,
          P@p,
          P@x,
          Y,
          theta,
          n_threads,
          gr_eps,
          eps
        )
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

bh_sparse_init <- function(
  perplexity,
  inp_kernel,
  symmetrize,
  normalize,
  n_threads,
  nn,
  row_normalize = TRUE,
  kernels = c("nngaussian", "knn")
) {
  inp_kernel <- match.arg(tolower(inp_kernel), kernels)
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
  function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
    cost <- sne_init(
      cost,
      X,
      perplexity = perplexity,
      kernel = inp_kernel,
      symmetrize = symmetrize,
      normalize = normalize,
      row_normalize = row_normalize,
      verbose = verbose,
      ret_extra = ret_extra,
      n_threads = n_threads,
      nn = nn,
      sparse = TRUE
    )
    cost$P <- Matrix::drop0(methods::as(cost$P, "generalMatrix"))
    cost
  }
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
  init <- bh_sparse_init(
    perplexity,
    inp_kernel,
    symmetrize,
    normalize,
    n_threads,
    nn
  )
  list(
    ndim = 2,
    init = function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
      cost <- init(cost, X, max_iter, verbose, ret_extra)
      cost$invN <- if (normalize) 1 else 1 / sum(cost$P)
      cost
    },
    cache_input = function(cost) {
      P <- cost$P
      if (!normalize) {
        # Keep the input-only term at the original affinity scale.
        if (is.null(cost$constant)) {
          cost$constant <- cost$invN * sum(P@x * log(P@x)) - 1 - log(lambda)
        }
        return(cost)
      }
      # Dense HSSNE floors the input-only terms. Count implicit
      # zeros analytically so that flooring never fills the sparse matrix.
      n <- nrow(P)
      values <- pmax(P@x, eps)
      missing <- n * (n - 1.0) - length(values)
      constant <- sum(values * log(values))
      if (eps > 0) {
        constant <- constant + missing * eps * log(eps)
      }
      cost$constant <- constant
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
