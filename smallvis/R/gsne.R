# global-SNE
gsne <- function(perplexity,
                 lambda = 1,
                 inp_kernel = "gaussian",
                 symmetrize = "symmetric",
                 eps = .Machine$double.eps,
                 use_cpp = FALSE,
                 n_threads = 0) {
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

      # Phat isn't affected by early exaggeration, so we cache it once only
      if (methods::is(X, "dist")) {
        Phat <- as.matrix(X)
      } else {
        Phat <- calc_d2(X, use_cpp = use_cpp, n_threads = n_threads)
      }

      Phat <- Phat + 1
      diag(Phat) <- 0

      Phat <- Phat / sum(Phat)
      cost$Phat <- Phat
      cost$phlogph <- colSums(Phat * logm(Phat, eps))

      cost
    },
    cache_input = function(cost) {
      eps <- cost$eps
      P <- cost$P
      P[P < eps] <- eps
      cost$plogp <- colSums(P * logm(P, eps))

      cost$plamphat <- P - lambda * cost$Phat

      cost
    },
    pfn = function(cost, Y) {
      cost <- cost_update(cost, Y)

      eps <- cost$eps
      kl <- cost$plogp - colSums(cost$P * logm(cost$W / cost$Z, eps))
      klhat <- cost$phlogph - colSums(cost$Phat * logm(cost$What / cost$Zhat, eps))
      cost$pcost <- kl + lambda * klhat
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)

      qlampqhat <- (cost$W / cost$Z) - lambda * (cost$What / cost$Zhat)
      cost$G <- k2g(Y, 4 * cost$W * (cost$plamphat - qlampqhat))
      cost
    },
    update = function(cost, Y) {
      W <- calc_d2(Y, use_cpp = use_cpp, n_threads = n_threads)
      What <- 1 + W
      W <- 1 / What
      diag(W) <- 0
      cost$Z <- sum(W)
      cost$W <- W

      diag(What) <- 0
      cost$What <- What
      cost$Zhat <- sum(What)

      cost
    }
  )
}
