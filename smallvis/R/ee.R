# Elastic Embedding -------------------------------------------------------

# Carreira-Perpinán, M. A. (2010, June).
# The Elastic Embedding Algorithm for Dimensionality Reduction.
# In \emph{Proceedings of the 27th International Conference on Machine Learning (ICML-10)} (pp. 167-174).
# http://faculty.ucmerced.edu/mcarreira-perpinan/papers/icml10.pdf (PDF)
# lambda control the strength of repulsive vs attractive forces
# if neg_weights is true, the repulsive contribution is weighted based on the
# squared input distances. Otherwise, no weighting is applied.
ee <- function(perplexity,
               lambda = 100,
               neg_weights = TRUE,
               inp_kernel = "gaussian",
               symmetrize = "symmetric",
               eps = .Machine$double.eps,
               use_cpp = FALSE,
               n_threads = 0) {
  list(
    init = function(cost,
                    X,
                    max_iter,
                    verbose = FALSE,
                    ret_extra = c()) {
      symmetrize <- match.arg(tolower(symmetrize), true_symmetrize_options())
      if (neg_weights) {
        if (methods::is(X, "dist")) {
          R <- X
        } else {
          R <- calc_d(X, use_cpp = use_cpp, n_threads = n_threads)
        }
        cost$Vn <- R / sum(R)
      } else {
        cost$Vn <- 1
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
      cost$G <- k2g(Y, 4 * (Vp - lambda * Vn * cost$W))
      cost
    },
    export = function(cost, val) {
      res <- NULL
      if (!is.null(cost[[val]])) {
        res <- cost[[val]]
      } else if (!is.null(cost[[toupper(val)]])) {
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


# t-Distributed Elastic Embedding
# EE-like cost function in terms of I-Divergence
# Scaled to give a gradient similar in form to t-SNE
tee <- function(perplexity,
                inp_kernel = "gaussian",
                symmetrize = "symmetric",
                lambda = 0.01,
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
      ret_extra <- unique(c(ret_extra, "V", "dint"))
      cost <- sne_init(
        cost,
        X,
        perplexity = perplexity,
        kernel = inp_kernel,
        symmetrize = symmetrize,
        normalize = FALSE,
        verbose = verbose,
        ret_extra = ret_extra,
        n_threads = n_threads,
        use_cpp = use_cpp
      )
      V <- cost$P
      cost$eps <- eps
      cost$invN <- 1 / sum(V)
      cost$gradconst <- 4 * cost$invN
      cost$lambda <- lambda

      V[V < eps] <- eps
      cost$constV <- cost$invN * (colSums(V * logm(V, eps)) - lambda * colSums(V))
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)
      cost$G <- k2g(Y, cost$gradconst * cost$W * (cost$P - cost$W * cost$lambda))
      cost
    },
    update = function(cost, Y) {
      W <- calc_d2(Y, use_cpp = use_cpp, n_threads = n_threads)
      W <- 1 / (1 + W)
      diag(W) <- 0
      cost$W <- W
      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)

      V <- cost$P
      W <- cost$W
      eps <- cost$eps

      cost$pcost <- cost$constV +
        cost$invN * (cost$lambda * colSums(W) - colSums(V * logm(W, eps)))
      cost
    },
    exaggerate = function(cost, exaggeration_factor) {
      cost$V <- cost$V * exaggeration_factor
      cost
    }
  )
}
