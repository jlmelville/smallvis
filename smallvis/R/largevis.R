# LargeVis ----------------------------------------------------------------

largevis <- function(perplexity,
                     inp_kernel = "gaussian",
                     symmetrize = "symmetric",
                     gamma = 1,
                     gr_eps = 0.1,
                     normalize = TRUE,
                     eps = 1e-9,
                     row_weight = NULL,
                     use_cpp = FALSE,
                     n_threads = 0) {
  if (!is.null(row_weight)) {
    row_normalize <- row_weight
  } else {
    row_normalize <- TRUE
  }
  lreplace(
    tsne(perplexity, use_cpp = use_cpp, n_threads = n_threads),
    init = function(cost,
                    X,
                    max_iter,
                    verbose = FALSE,
                    ret_extra = c()) {
      symmetrize <- match.arg(tolower(symmetrize), true_symmetrize_options())
      cost <- sne_init(
        cost,
        X = X,
        perplexity = perplexity,
        symmetrize = symmetrize,
        kernel = inp_kernel,
        normalize = normalize,
        verbose = verbose,
        row_normalize = row_normalize,
        ret_extra = ret_extra,
        n_threads = n_threads,
        use_cpp = use_cpp
      )
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
      cost$G <- k2g(Y, 4 * W * (cost$P - ((gamma * W) / (
        1 + cost$greps1 * W
      ))))
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
