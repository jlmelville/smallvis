bh_ssne <- function(
  perplexity,
  inp_kernel = "nngaussian",
  symmetrize = "symmetric",
  n_threads = 0,
  theta = 0.5,
  nn = "approximate"
) {
  bh_sne(
    perplexity,
    inp_kernel,
    symmetrize,
    eps = 0,
    n_threads = n_threads,
    theta = theta,
    nn = nn,
    alpha = 0,
    normalize = TRUE,
    lambda = 0
  )
}

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

bhdhssne <- function(
  perplexity,
  alpha = 0.5,
  inp_kernel = "nngaussian",
  symmetrize = "symmetric",
  eps = .Machine$double.eps,
  n_threads = 0,
  theta = 0.5,
  nn = "approximate"
) {
  lreplace(
    bhhssne(
      perplexity,
      alpha,
      inp_kernel,
      symmetrize,
      eps,
      n_threads,
      theta,
      nn
    ),
    gradient_only = TRUE,
    epoch = dhssne_epoch
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

bhjssne <- function(
  perplexity,
  inp_kernel = "nngaussian",
  symmetrize = "symmetric",
  eps = .Machine$double.eps,
  n_threads = 0,
  theta = 0.5,
  nn = "approximate"
) {
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("Barnes-Hut Jensen-Shannon SNE requires a finite positive eps")
  }
  bh_cost(
    init = bh_sparse_init(
      perplexity,
      inp_kernel,
      symmetrize,
      normalize = TRUE,
      n_threads = n_threads,
      nn = nn
    ),
    cache_input = function(cost) {
      P <- cost$P
      p <- pmax(P@x, eps)
      # Match the dense entropy floor on implicit off-diagonal zeros.
      off_graph <- nrow(P) * (nrow(P) - 1) - length(p)
      cost$plogp <- sum(p * log(p)) + off_graph * eps * log(eps)
      cost
    },
    pfn = function(cost, Y) {
      P <- cost$P
      score <- bh_jssne_score_cpp(P@i, P@p, P@x, Y, theta, n_threads)
      cost$pcost <- 0.5 * (cost$plogp + score)
      cost
    },
    gr = function(cost, Y) {
      P <- cost$P
      cost$G <- bh_jssne_gradient_cpp(P@i, P@p, P@x, Y, theta, n_threads)
      cost
    }
  )
}

bhhlsne <- function(
  perplexity,
  inp_kernel = "nngaussian",
  n_threads = 0,
  theta = 0.5,
  nn = "approximate"
) {
  bh_cost(
    init = bh_sparse_init(
      perplexity,
      inp_kernel,
      symmetrize = "symmetric",
      normalize = TRUE,
      n_threads = n_threads,
      nn = nn
    ),
    cache_input = function(cost) {
      p <- cost$P@x
      cost$sqrt_p <- sqrt(p)
      cost$psum <- sum(p)
      cost
    },
    pfn = function(cost, Y) {
      P <- cost$P
      similarity <- bh_hlsne_similarity_cpp(
        P@i,
        P@p,
        cost$sqrt_p,
        Y,
        theta,
        n_threads
      )
      cost$pcost <- cost$psum + 1 - 2 * similarity
      cost
    },
    gr = function(cost, Y) {
      P <- cost$P
      cost$G <- bh_hlsne_gradient_cpp(
        P@i,
        P@p,
        cost$sqrt_p,
        Y,
        theta,
        n_threads
      )
      cost
    }
  )
}

bhchsne <- function(
  perplexity,
  inp_kernel = "nngaussian",
  symmetrize = "symmetric",
  n_threads = 0,
  theta = 0.5,
  nn = "approximate"
) {
  bh_cost(
    init = bh_sparse_init(
      perplexity,
      inp_kernel,
      symmetrize,
      normalize = TRUE,
      n_threads = n_threads,
      nn = nn
    ),
    cache_input = function(cost) {
      p <- cost$P@x
      cost$p_squared <- p * p
      cost$constant <- 1 - 2 * sum(p)
      cost
    },
    pfn = function(cost, Y) {
      P <- cost$P
      cost$pcost <- cost$constant +
        bh_chsne_score_cpp(P@i, P@p, cost$p_squared, Y, theta, n_threads)
      cost
    },
    gr = function(cost, Y) {
      P <- cost$P
      cost$G <- bh_chsne_gradient_cpp(
        P@i,
        P@p,
        cost$p_squared,
        Y,
        theta,
        n_threads
      )
      cost
    }
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
  bh_cost(
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
    }
  )
}

bhumap <- function(
  perplexity,
  inp_kernel = "skd",
  symmetrize = "umap",
  spread = 1,
  min_dist = 0.001,
  gr_eps = 0.1,
  eps = 1e-9,
  row_weight = NULL,
  n_threads = 0,
  theta = 0.5,
  nn = "approximate"
) {
  bh_bernoulli(
    perplexity,
    inp_kernel,
    symmetrize,
    gr_eps,
    eps,
    row_weight,
    n_threads,
    theta,
    nn,
    spread = spread,
    min_dist = min_dist,
    student = FALSE
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
  bh_bernoulli(
    perplexity,
    inp_kernel,
    symmetrize,
    gr_eps,
    eps,
    row_weight,
    n_threads,
    theta,
    nn,
    student = TRUE
  )
}

bh_bernoulli <- function(
  perplexity,
  inp_kernel,
  symmetrize,
  gr_eps,
  eps,
  row_weight,
  n_threads,
  theta,
  nn,
  student,
  spread = 1,
  min_dist = 0.001
) {
  label <- if (student) "t-UMAP" else "UMAP"
  if (
    !is.numeric(gr_eps) ||
      length(gr_eps) != 1L ||
      !is.finite(gr_eps) ||
      gr_eps < 0
  ) {
    stop("Barnes-Hut ", label, " requires a finite nonnegative gr_eps")
  }
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("Barnes-Hut ", label, " requires a finite positive eps")
  }
  init <- bh_sparse_init(
    perplexity,
    inp_kernel,
    symmetrize,
    normalize = FALSE,
    n_threads = n_threads,
    nn = nn,
    row_normalize = if (is.null(row_weight)) FALSE else row_weight,
    kernels = if (student) c("nngaussian", "knn", "skd") else "skd"
  )
  bh_cost(
    gradient_only = gr_eps != 0,
    init = if (student) {
      init
    } else {
      function(cost, X, max_iter, verbose = FALSE, ret_extra = c()) {
        cost <- init(cost, X, max_iter, verbose, ret_extra)
        init_ab(cost, spread = spread, min_dist = min_dist, verbose = verbose)
      }
    },
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
      score <- if (student) {
        bh_tumap_score_cpp(
          P@i,
          P@p,
          P@x,
          Y,
          theta,
          n_threads,
          eps
        )
      } else {
        bh_umap_score_cpp(
          P@i,
          P@p,
          P@x,
          Y,
          theta,
          n_threads,
          eps,
          cost$a,
          cost$b
        )
      }
      cost$pcost <- cost$Cp - score
      cost
    },
    gr = function(cost, Y) {
      P <- cost$P
      G <- if (student) {
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
      } else {
        bh_umap_gradient_cpp(
          P@i,
          P@p,
          P@x,
          Y,
          theta,
          n_threads,
          gr_eps,
          eps,
          cost$a,
          cost$b
        )
      }
      cost$G <- 4 * G
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
  bh_cost(
    alpha = alpha,
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
        cost$alpha,
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
        cost$alpha,
        normalize,
        lambda
      )
      cost$G <- (4 * cost$invN) * G
      cost
    }
  )
}

bhrklsne <- function(
  perplexity,
  inp_kernel = "nngaussian",
  symmetrize = "symmetric",
  eps = .Machine$double.eps,
  n_threads = 0,
  theta = 0.5,
  nn = "approximate"
) {
  bh_mixed_kl(
    perplexity,
    inp_kernel,
    symmetrize,
    eps,
    n_threads,
    theta,
    nn,
    lambda = 0,
    gaussian = FALSE
  )
}

bhsnerv <- function(
  perplexity,
  lambda = 0.9,
  inp_kernel = "nngaussian",
  symmetrize = "symmetric",
  eps = .Machine$double.eps,
  n_threads = 0,
  theta = 0.5,
  nn = "approximate"
) {
  bh_mixed_kl(
    perplexity,
    inp_kernel,
    symmetrize,
    eps,
    n_threads,
    theta,
    nn,
    lambda,
    gaussian = TRUE
  )
}

bh_mixed_kl <- function(
  perplexity,
  inp_kernel,
  symmetrize,
  eps,
  n_threads,
  theta,
  nn,
  lambda,
  gaussian
) {
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("Barnes-Hut KL requires a finite positive eps")
  }
  if (
    !is.numeric(lambda) ||
      length(lambda) != 1L ||
      !is.finite(lambda) ||
      lambda < 0 ||
      lambda > 1
  ) {
    stop("Barnes-Hut SNeRV requires lambda between 0 and 1")
  }
  evaluate <- function(cost, Y, gradient) {
    P <- cost$P
    bh_mixed_kl_cpp(
      P@i,
      P@p,
      cost$log_ratio,
      cost$p_excess,
      Y,
      theta,
      n_threads,
      eps,
      lambda,
      cost$p_mass,
      cost$plogp,
      gaussian,
      gradient
    )
  }
  bh_cost(
    init = bh_sparse_init(
      perplexity,
      inp_kernel,
      symmetrize,
      normalize = TRUE,
      n_threads = n_threads,
      nn = nn
    ),
    cache_input = function(cost) {
      P <- cost$P
      p <- pmax(P@x, eps)
      log_p <- log(p)
      n <- nrow(P)
      off_graph <- n * (n - 1.0) - length(p)
      cost$log_ratio <- log_p - log(eps)
      cost$p_excess <- p - eps
      cost$p_mass <- sum(p) + off_graph * eps
      cost$plogp <- sum(p * log_p) + off_graph * eps * log(eps)
      cost
    },
    pfn = function(cost, Y) {
      cost$pcost <- evaluate(cost, Y, gradient = FALSE)
      cost
    },
    gr = function(cost, Y) {
      cost$G <- evaluate(cost, Y, gradient = TRUE)
      cost
    }
  )
}

bh_cost <- function(...) {
  list(
    ndim = 2,
    ...,
    export = cost_export,
    update = function(cost, Y) cost,
    exaggerate = function(cost, exaggeration_factor) {
      cost$P <- cost$P * exaggeration_factor
      cost
    }
  )
}
