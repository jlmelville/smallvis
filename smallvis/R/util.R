# df_sub <- function(df, patterns) {
#   if (length(patterns) == 0) {
#     return(df)
#   }
#   if (length(patterns) == 1) {
#     df[, grepl(patterns, names(df))]
#   }
#   else {
#     l <- list()
#     for (pattern in patterns) {
#       l[[pattern]] <- df[, grepl(pattern, names(df))]
#     }
#     do.call("cbind", l)
#   }
# }
#
#
# # # Clamp Numerical Values
# # #
# # # Values are truncated so that they lie within (\code{min_val, max_val}). In
# # # embedding this is used to prevent individual probabilities values getting
# # # too small and causing underflow or some other horrible explosion.
# # #
# # # @param x Matrix.
# # # @param min_val Minimum value allowed for any element in the matrix.
# # # @param max_val Maximum value allowed for any element in the matrix.
# # # @return Matrix with the clamped values.
# # clamp <- function(x, min_val = .Machine$double.eps, max_val = NULL) {
# #   x[x < min_val] <- min_val
# #   if (!is.null(max_val)) {
# #     x[x > max_val] <- max_val
# #   }
# #   x
# # }
#
# tsne_cost <- function(P, Y, eps = .Machine$double.eps) {
#   Q <- dist2(Y)
#   Q <- 1 / (1 + Q)
#   diag(Q) <- 0
#   Q <- Q / sum(Q)
#   Q[Q < eps] <- eps
#
#   sum(P * log((P + eps) / (Q + eps)))
# }
#
gradient_fd <- function(Y, cost_fn, diff = .Machine$double.eps ^ (1 / 3)) {
  cost <- function(Y) {
    cost_fn <- cost_fn$gr(cost_fn, Y)
    cost_fn <- cost_fn$pfn(cost_fn, Y)
    sum(cost_fn$pcost)
  }

  nr <- nrow(Y)
  nc <- ncol(Y)

  Gfd <- matrix(0, nrow = nr, ncol = nc)
  for (i in 1:nr) {
    for (j in 1:nc) {
      Yij_old <- Y[i, j]

      Y[i, j] <- Yij_old + diff
      cost_fwd <- cost(Y)

      Y[i, j] <- Yij_old - diff
      cost_back <- cost(Y)

      fd <- (cost_fwd - cost_back) / (2 * diff)
      Gfd[i, j] <- fd

      Y[i, j] <- Yij_old
    }
  }

  Gfd
}
#
# dot <- function(X, Y = X) {
#   sum(X * Y)
# }
#
# smooth_knn_distances_ref <- function(distances, k, n_iter = 128,
#                                  smooth_k_tolerance = 1e-5,
#                                  min_k_dist_scale = 1e-3) {
#   k <- min(k, nrow(distances) - 1)
#   target <- log2(k)
#   rho <- rep(0, nrow(distances))
#   result <- rep(0, nrow(distances))
#
#   for (i in 1:nrow(distances)) {
#     lo <- 0.0
#     hi <- Inf
#     mid <- 1.0
#     distances[i, ] <- sort(distances[i, ])
#     # TODO: This is very inefficient, but will do for now. FIXME
#     ith_distances <- distances[i, ]
#     non_zero_dists <- ith_distances[ith_distances > 0.0]
#
#     if (length(non_zero_dists) > 0) {
#       rho[i] <- min(non_zero_dists)
#     }
#     else {
#       rho[i] <- 0.0
#     }
#     for (n in 1:n_iter) {
#       psum <- 0.0
#       for (j in 2:(k + 1)) {
#         psum <- psum + exp(-((distances[i, j] - rho[i]) / mid))
#       }
#
#       val <- psum
#
#       if (abs(val - target) < smooth_k_tolerance) {
#         break
#       }
#
#       if (val > target) {
#         hi <- mid
#         mid <- (lo + hi) / 2.0
#       }
#       else {
#         lo <- mid
#         if (is.infinite(hi)) {
#           mid <- mid * 2
#         }
#         else {
#           mid <- (lo + hi) / 2.0
#         }
#       }
#
#       if (i == 143) {
#         message(n_iter, " ", mid, " ", val, " ", abs(val - target), " rho = ", rho[i])
#         browser()
#       }
#     }
#     result[i] <- mid
#
#     # TODO: This is very inefficient, but will do for now. FIXME
#     if (rho[i] > 0.0) {
#       if (result[i] < min_k_dist_scale * mean(ith_distances)) {
#         result[i] <- min_k_dist_scale * mean(ith_distances)
#       }
#     }
#     else {
#       if (result[i] < min_k_dist_scale * mean(distances)) {
#         result[i] <- min_k_dist_scale * mean(distances)
#       }
#     }
#   }
#   list(
#     sigmas = result,
#     rhos = rho
#   )
# }
#
# plotumapfun <- function(spread = 1, min_dist = 0.001) {
#   xv <- seq(from = 0, to = spread * 3, length.out = 300)
#   yv <- rep(0, length(xv))
#   yv[xv < min_dist] <- 1
#   yv[xv >= min_dist] <- exp(-(xv[xv >= min_dist] - min_dist) / spread)
#   graphics::plot(xv, yv, type = "l")
# # stats::nls(yv ~ 1 / (1 + a * xv ^ (2 * b)),
# #            start = list(a = 1, b = 1))$m$getPars()
# }
#
embed_img <- function(df, res, title = "", perplexity = 40, equal_axes = FALSE,
                      pc_axes = FALSE,
                      ref = NULL, xflip = FALSE, yflip = FALSE, rot = NULL,
                      title_cb = mnp_cb, sub = "", cex = 1) {
  if (methods::is(res, "matrix")) {
    res <- list(Y = res)
  }
  if (is.null(res$DX)) {
    indexes <- which(vapply(df, is.numeric, logical(1)))
    X <- as.matrix(df[, indexes])
    res$DX <- sqrt(safe_dist2(X))
  }

  if (is.null(res$perplexity)) {
    res$perplexity <- perplexity
  }

  if (!is.null(res$alpha)) {
    title <- paste0(title, "(", formatC(res$alpha), ")")
  }

  if (!is.null(title_cb)) {
    title <- paste0(title, title_cb(res))
  }

  if (!is.null(ref)) {
    res$Y <- kabsch(ref$Y, res$Y)
  }
  if (xflip) {
    res$Y <- flipX(res$Y)
  }
  if (yflip) {
    res$Y <- flipY(res$Y)
  }
  if (!is.null(rot)) {
    res$Y <- rotate2D(res$Y, rot)
  }
  vizier::embed_plot(res$Y, df, title = title, equal_axes = equal_axes,
                     pc_axes = pc_axes, sub = sub, cex = cex)
}

mnp_cb <- function(res) {
  if (is.null(res$DY)) {
    res$DY <- sqrt(safe_dist2(res$Y))
  }
  av_pres <- quadra::nbr_pres(res$DX, res$DY, res$perplexity)
  paste0(" mnp@", formatC(res$perplexity), " = ", formatC(mean(av_pres)))
}

# Kabsch Algorithm
#
# Aligns two sets of points via rotations and translations.
#
# Given two sets of points, with one specified as the reference set,
# the other set will be rotated so that the RMSD between the two is minimized.
# The format of the matrix is that there should be one row for each of
# n observations, and the number of columns, d, specifies the dimensionality
# of the points. The point sets must be of equal size and with the same
# ordering, i.e. point one of the second matrix is mapped to point one of
# the reference matrix, point two of the second matrix is mapped to point two
# of the reference matrix, and so on.
#
# param pm n x d matrix of reference points.
# param qm n x d matrix of points to align to to \code{pm}
# return Matrix \code{qm} rotated and translated so that the ith point
#  is aligned to the ith point of \code{pm} in the least-squares sense.
# references
# \url{https://en.wikipedia.org/wiki/Kabsch_algorithm}
kabsch <- function(pm, qm) {
  pm_dims <- dim(pm)
  if (!all(dim(qm) == pm_dims)) {
    stop(call. = TRUE, "Point sets must have the same dimensions")
  }
  # The rotation matrix will have (ncol - 1) leading ones in the diagonal
  diag_ones <- rep(1, pm_dims[2] - 1)

  # center the points
  pm <- scale(pm, center = TRUE, scale = FALSE)
  qm <- scale(qm, center = TRUE, scale = FALSE)

  am <- crossprod(pm, qm)

  svd_res <- svd(am)
  # use the sign of the determinant to ensure a right-hand coordinate system
  d <- determinant(tcrossprod(svd_res$v, svd_res$u))$sign
  dm <- diag(c(diag_ones, d))

  # rotation matrix
  um <- svd_res$v %*% tcrossprod(dm, svd_res$u)

  # Rotate and then translate to the original centroid location of pm
  sweep(t(tcrossprod(um, qm)), 2, -attr(pm, "scaled:center"))
}

flipY <- function(X) {
  X[, 2] <- -X[, 2]
  X
}

flipX <- function(X) {
  X[, 1] <- -X[, 1]
  X
}

rotate2D <- function(X, deg = 90) {
  theta <- deg * (pi / 180)
  R <- matrix(c(
    cos(theta), -sin(theta),
    sin(theta), cos(theta)),
    byrow = TRUE, nrow = 2
  )

  t(R %*% t(X))
}
#
# umapfn <- function(P, Y, eps = .Machine$double.eps,
#                    spread = 1, min_dist = 0.001) {
#   ab_params <- find_ab_params(spread = spread, min_dist = min_dist)
#   a <- ab_params[1]
#   b <- ab_params[2]
#   W <- dist2(Y)
#   W <- 1 / (1 + a * W ^ b)
#   sum(-P * log(W + eps) - (1 - P) * log1p(-W + eps))
# }
#
# tsnefn <- function(P, Y, eps = .Machine$double.eps) {
#   W <- dist2(Y)
#   W <- 1 / (1 + W)
#   W <- W / sum(W)
#   sum(P * log((P + eps)/(W + eps)))
# }
#
# lvfn <- function(P, Y, eps = .Machine$double.eps, gamma = 7) {
#   W <- dist2(Y)
#   W <- 1 / (1 + W)
#   sum(-P * log(W + eps) - gamma * log1p(-W + eps))
# }
#
# tumapfn <- function(P, Y, eps = .Machine$double.eps) {
#   W <- dist2(Y)
#   W <- 1 / (1 + W)
#   sum(-P * log(W + eps) - (1 - P) * log1p(-W + eps))
# }
#
# tumaprepfn <- function(P, Y, eps = .Machine$double.eps) {
#   D2 <- dist2(Y)
#   W <- 1 / (1 + D2)
#   sum((1 - P) * log((1 - W) + eps))
# }
#
# rmsd <- function(X, Y) {
#   sqrt(sum((X - Y) ^ 2) / nrow(X))
# }
#
# specd <- function(P, G) {
#   # do this once
#   L <- diag(colSums(P)) - P
#   mu <- min(L[L > 0]) * 1e-10
#   H <- 4 * (L + mu)
#   R <- chol(H)
#
#   # solve once per G
#   backsolve(R, backsolve(R, -G, transpose = TRUE))
# }
#
# remove_nulls <- function(l) {
#   l[!vapply(l, is.null, logical(1))]
# }
#
#
#
#
# # else if (opt == "spectral") {
# #   if (iter == 1) {
# #     # L
# #     R <- diag(colSums(P)) - P
# #     mu <- min(R[R > 0]) * 1e-10
# #     # H
# #     R <- 4 * (R + mu)
# #     # cholesky decomposition of H
# #     R <- chol(R)
# #   }
# #
# #   uY <- backsolve(R, backsolve(R, -G, transpose = TRUE))
# # }
#
# # // [[Rcpp::export]]
# # NumericMatrix dist2c(const NumericMatrix& x) {
# #   const int N = x.nrow();
# #   const int K = x.ncol();
# #   NumericMatrix d2(N, N);
# #
# #   double ddsum = 0.0;
# #   double dd = 0.0;
# #   for (int i = 0; i < N; i++) {
# #     for (int j = i + 1; j < N; j++) {
# #       ddsum = 0.0;
# #       for (int k = 0; k < K; k++) {
# #         dd = x(i, k) - x(j, k);
# #         ddsum += dd * dd;
# #       }
# #       d2(i, j) = ddsum;
# #       d2(j, i) = ddsum;
# #     }
# #   }
# #
# #   return d2;
# # }
#
#
# # else if (method == "ntumap") {
# #   Q <- W / sum(W)
# #   C <- (P - Q) / (1 - Q)
# #   G <- 4 * W * ((P - Q) / (1 - Q) - sum(C) * Q)
# # }
#
#
#
# # Floyd's algorithm for finding the shortest paths between two nodes,
# # starting from the distance matrix, with nearest neighbor distances marked
# # and everything else Inf
# # Unlike the Rfast version which works directly on the
# # slow_floyd <- function(D, verbose = FALSE) {
# #   n <- nrow(D)
# #
# #   n1s <- rep(1, each = n)
# #   for (i in 1:n) {
# #     di <- t(D[, i])
# #     di <- t(di[n1s, ])
# #     tdi <- t(D[i, ])
# #     tdi <- tdi[n1s, ]
# #     D <- pmin(D, di + tdi)
# #     if (verbose && i %% floor(n / 10) == 0) {
# #       message("Calculated shortest paths for ", i, " / ", n)
# #     }
# #   }
# #   D
# # }
#
#

loadsnedata <- function() {
  load(file = "/dev/R/datasets/snedata.Rda", envir = .GlobalEnv)
  load(file = "/dev/R/datasets/sr3k.Rda", envir = .GlobalEnv)
}

repm <- function(X, n) {
  nr <- nrow(X)
  nc <- ncol(X)
  nnr <- nr * n
  nnc <- nc * n

  m <- matrix(0, nrow = nnr, ncol = nnc)
  for (i in 1:n) {
    m[(((i - 1) * nr) + 1):(nr * i), (((i - 1) * nc) + 1):(nc * i)] <- X
  }
  m
}

pca_plot <- function(df, title = NULL) {
  indexes <- which(vapply(df, is.numeric, logical(1)))
  X <- as.matrix(df[, indexes])
  vizier::embed_plot(pca_scores(X, ncol = 2), df, title = title)
}

dot <- function(a, b = a) {
  sum(a * a)
}

rmsy <- function(X, Y) {
  dx <- sqrt(safe_dist2(X))
  dy <- sqrt(safe_dist2(Y))

  sqrt(0.5 * sum((dx - dy) ^ 2)) / nrow(X)

}

# convert data frame to matrix using numeric columns
x2m <- function(X) {
  as.matrix(X[, which(vapply(X, is.numeric, logical(1)))])
}


# More expensive but more generic intrinsic dimensionality calculation
# where only a vector of exponential affinities, Wi, need to be available
intrinsic_dimensionality <- function(Wi, eps = .Machine$double.eps) {
  iZ <- 1 / sum(Wi)
  lw <- log(Wi + eps)

  wlwlw <- sum(Wi * lw * lw)

  wlw2 <- sum(Wi * lw)
  wlw2 <- wlw2 * wlw2

  2 * iZ * (wlwlw - iZ * wlw2)
}

intdim_exp <- function(W, eps = .Machine$double.eps) {
  n <- nrow(W)
  id <- rep(0, n)
  for (i in 1:n) {
    id[i] <- intrinsic_dimensionality(W[i, ], eps = eps)
  }

  id
}

# intdim_Wt <- function(W, eps = .Machine$double.eps) {
#   n <- nrow(W)
#   id <- rep(0, n)
#
#   for (i in 1:n) {
#     id[i] <- intdim_tdist(W[i, ], eps = eps)
#   }
#
#   id
# }
#
# intdim_tdist <- function(Wi, eps = .Machine$double.eps) {
#   Z <- sum(Wi)
#   iz2 <- 1 / (Z * Z)
#
#   lw <- log(Wi + eps)
#   wlw <- sum(Wi * lw)
#
#   inner <- Z * lw - wlw
#
#   2 * iz2 * sum((Wi - 1) * Wi * inner)
# }
#
# intdim_k <- function(X, k) {
#   knn <- FNN::get.knn(X, k = k)
#
#   dx <- log(k) - log(k - 1)
#   dy <- log(knn$nn.dist[, k]) - log(knn$nn.dist[, k - 1])
#
#   summarize(dx / dy)
# }
#
# intdim_knn <- function(X, k = nrow(X) - 1, base = exp(1)) {
#   knn <- FNN::get.knn(X, k)
#   dlrad <- log(knn$nn.dist[, k] / knn$nn.dist[, k - 1], base = base)
#   dlk <- log(k / (k - 1), base = base)
#
# browser()
#   stats::lm(log(1:k) ~ 0 + log(knn$nn.dist[1, 1:k]))
#   #
#   #
#   # dlk / dlrad
# }
#
#
# # Facco, E., d’Errico, M., Rodriguez, A., & Laio, A. (2017).
# # Estimating the intrinsic dimension of datasets by a minimal neighborhood
# # information.
# # Scientific reports, 7(1), 12140.
# twonn <- function(X, eps = .Machine$double.eps) {
#   knn <- FNN::get.knn(X, 2)
#   mu <- knn$nn.dist[, 2] / (knn$nn.dist[, 1] + eps)
#
#   mu <- sort(mu)
#
#   Femp <- (0:(length(mu) - 1)) / length(mu)
#   p90 <- ceiling(length(mu) * 0.9)
#   mu <- mu[1:p90]
#   Femp <- Femp[1:p90]
#
#   plot(log(mu), -log1p(-Femp))
#   lm(-log1p(-Femp) ~ log(mu))$coefficients[2]
# }
#
# twonnp <- function(X, eps = .Machine$double.eps) {
#   knn <- FNN::get.knn(X, 2)
#   mu <- knn$nn.dist[, 2] / (knn$nn.dist[, 1] + eps)
#
#   mu <- sort(mu)
#
#   Femp <- (0:(length(mu) - 1)) / length(mu)
#   p90 <- ceiling(length(mu) * 0.9)
#   mu <- mu[1:p90]
#   Femp <- Femp[1:p90]
#
#   plot(log(mu), -log1p(-Femp))
#   lm(-log1p(-Femp) ~ log(mu))$coefficients[2]
# }
#
# intdim_knnc <- function(X, k = nrow(X) - 1, base = exp(1), eps = .Machine$double.eps) {
#   knn <- FNN::get.knn(X, k + 1)
#   dx <- log(knn$nn.dist[, k + 1] / knn$nn.dist[, k - 1], base = base)
#   dy <- log((k + 1) / (k - 1), base = base)
#   dy / dx
# }
#
# intdim_knn_k <- function(knn, k) {
#   dy <- log(k) - log(k - 1)
#   dx <- log(knn$nn.dist[, k]) - log(knn$nn.dist[, k - 1])
#
#   dy / dx
# }
#
# radk <- function(X, base = exp(1)) {
#   n <- nrow(X)
#   D <- sqrt(safe_dist2(X))
#   D <- t(apply(D, 1, sort))
#   meanD <- colMeans(D)
#   mediD <- apply(D, 2, median)
#   statD <- mediD
#   grads <- rep(0, n - 1)
#   for (i in 1:(n - 2)) {
#     dy <- log(i + 1, base = base) - log(i, base = base)
#     dx <- log(statD[i + 2], base = base) - log(statD[i + 1], base = base)
#     grads[i] <- dy / dx
#   }
#   statD <- statD[2:length(statD)]
#
#   plot(log(statD, base = base), log(1:(n - 1), base = base))
#   median(grads)
# }
#
# corr_int <- function(X) {
#   D <- sqrt(safe_dist2(X))
#   D <- sort(D[upper.tri(D)])
#   n <- length(D)
#   graphics::plot(log(D), log(1 / (1:n)))
# }
#
# int_dim_knn_av <- function(xfact, ks = c(2, 4, 8, 16, 32, 64, 128, 256, 512),
#                            nrep = 100) {
#   intdims <- matrix(nrow = nrep, ncol = length(ks))
#   for (i in 1:nrep) {
#     x <- x2m(xfact(n = 10000, dim = 1))
#     for (k in 1:length(ks)) {
#       intdims[i, k] <- stats::median(intdim_knnc(x, k = ks[k]))
#     }
#   }
#   colMeans(intdims)
# }
#
# ucube <- function(n = 1000, dim = 3, edge_length = 1) {
#   vec <- seq(from = 0, to = edge_length, length.out = ceiling(n ^ (1 / dim)))
#   lst <- lapply(numeric(dim), function(x) vec)
#   data.frame(expand.grid(lst))
# }
#
# jucubex <- function(n = 1000, dim = 3, edge_length = 1, factor = 1e-4) {
#   jitter(x2m(ucube(n = n, dim = dim, edge_length = edge_length)), factor = factor)
# }
#
#
# xtsne <- function(perplexity, inp_kernel = "gaussian") {
#   list(
#     init = function(cost, X, eps = .Machine$double.eps, verbose = FALSE,
#                     ret_extra = c()) {
#       ret_extra <- unique(c(ret_extra, "V"))
#       cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
#                        symmetrize = "symmetric", normalize = TRUE,
#                        verbose = verbose, ret_extra = ret_extra)
#       cost$eps <- eps
#       V <- cost$V
#       V <- V / rowSums(V)
#       V <- 0.5 * (V + t(V))
#       cost$V <- V
#       cost$sumV <- sum(V)
#       cost$isumV <- 1 / cost$sumV
#       cost
#     },
#     pfn = kl_cost,
#     gr = function(cost, Y) {
#       isumV <- cost$isumV
#       sumV <- cost$sumV
#       V <- cost$V
#       W <- dist2(Y)
#       W <- 1 / (1 + W)
#       diag(W) <- 0
#       invZ <- 1 / sum(W)
#       cost$invZ <- invZ
#       cost$W <- W
#       message(" rep weight = ", formatC(sumV * invZ))
#       cost$G <- k2g(Y, 4 * isumV * W * (V - W * sumV * invZ))
#       cost
#     },
#     export = function(cost, val) {
#       res <- NULL
#       if (!is.null(cost[[val]])) {
#         res <- cost[[val]]
#       }
#       else if (!is.null(cost[[toupper(val)]])) {
#         res <- cost[[toupper(val)]]
#       }
#       else {
#         switch(val,
#                q = {
#                  res <- cost$W * cost$invZ
#                })
#       }
#       res
#     }
#   )
# }

bsvbench <- function(...) {
  res <- list()
  varargs <- list(...)

  # res$iris <- do.call(smallvis, lreplace(X = iris, varargs))
  # res$s1k <- do.call(smallvis, lreplace(X = s1k, varargs))
  res$oli <- do.call(smallvis, lreplace(X = oli, varargs))
  res$frey <- do.call(smallvis, lreplace(X = frey, varargs))
  # res$coil20 <- do.call(smallvis, lreplace(X = coil20, varargs))
  res$mnist <- do.call(smallvis, lreplace(X = mnist6k, varargs))
  res$fashion <- do.call(smallvis, lreplace(X = fashion6k, varargs))

  res$varargs <- varargs

  res
}

bsvimg <- function(res, title = "") {
  # embed_img(iris, res$iris, title)
  # embed_img(s1k, res$s1k, title)
  embed_img(oli, res$oli, title)
  embed_img(frey, res$frey, title)
  # embed_img(coil20, res$coil20, title)
  embed_img(mnist6k, res$mnist, title)
  embed_img(fashion6k, res$fashion, title)
}

svbench <- function(img_only = FALSE, inc_swiss = FALSE, ...) {
  res <- list()
  varargs <- list(...)

  if (!img_only) {
    res$iris <- do.call(smallvis, lreplace(X = datasets::iris, varargs))
    res$s1k <- do.call(smallvis, lreplace(X = s1k, varargs))
  }
  res$oli <- do.call(smallvis, lreplace(X = oli, varargs))
  res$frey <- do.call(smallvis, lreplace(X = frey, varargs))
  res$coil20 <- do.call(smallvis, lreplace(X = coil20, varargs))
  res$mnist <- do.call(smallvis, lreplace(X = mnist6k, varargs))
  res$fashion <- do.call(smallvis, lreplace(X = fashion6k, varargs))
  if (inc_swiss) {
    res$sr3k <- do.call(smallvis, lreplace(X = sr3k, varargs))
  }

  res$varargs <- varargs

  res
}

svimg <- function(res, title = "", img_only = FALSE, inc_swiss = FALSE,
                  title_cb = mnp_cb, sub = "", perplexity = 40) {
  if (!img_only) {
    embed_img(datasets::iris, res$iris, title, title_cb = title_cb, sub = sub,
              perplexity = perplexity)
    embed_img(s1k, res$s1k, title, title_cb = title_cb, sub = sub,
              perplexity = perplexity)
  }
  embed_img(oli, res$oli, title, title_cb = title_cb, sub = sub,
            perplexity = perplexity)
  embed_img(frey, res$frey, title, title_cb = title_cb, sub = sub,
            perplexity = perplexity)
  embed_img(coil20, res$coil20, title, title_cb = title_cb, sub = sub,
            perplexity = perplexity)
  embed_img(mnist6k, res$mnist, title, title_cb = title_cb, sub = sub,
            perplexity = perplexity)
  embed_img(fashion6k, res$fashion, title, title_cb = title_cb, sub = sub,
            perplexity = perplexity)
  if (inc_swiss) {
    embed_img(sr3k, res$sr3k, title, title_cb = title_cb, sub = sub,
              perplexity = perplexity)
  }
}

clustbench <- function(...) {
  res <- list()
  varargs <- list(...)

  for (perp in c(2, 5, 30, 50, 100)) {
    res[[perp]] <- do.call(smallvis, lreplace(X = threec_50,
                                              perplexity = perp, varargs))
  }
  res$varargs <- varargs

  res
}

pm2tc <- function(X, digits = 4) {
 X <- signif(X, digits = digits)
 dim(X) <- NULL
 dput(X)
}

best_repn <- function(all_res) {
  nres <- length(all_res)
  best_res <- NULL
  best_cost <- Inf
  all_costs <- rep(0, nres)
  for (i in 1:nres) {
    res <- all_res[[i]]
    final_cost <- res$itercosts[length(res$itercosts)]
    if (final_cost < best_cost) {
      best_cost <- final_cost
      best_res <- res
    }
    names(final_cost) <- NULL
    all_costs[i] <- final_cost
  }
  best_res$all_costs <- all_costs
  best_res
}

# as.dist(combine_dy(iris_repn))
combine_dy <- function(all_res, stat_fun = stats::median) {
  res_dy <- NULL
  nres <- length(all_res)
  resi <- all_res[[1]]
  if (class(resi) == "matrix") {
    Y <- resi
  }
  else {
    Y <- resi$Y
  }
  nr <- nrow(Y)
  dijs <- rep(0, nres)

  res_dy <- matrix(0, nrow = nr, ncol = nr)
  for (i in 1:(nr - 1)) {
    for (j in (i + 1):nr) {
      for (k in 1:nres) {
        if (class(resi) == "matrix") {
          Y <- all_res[[k]]
        }
        else {
          Y <- all_res[[k]]$Y
        }
        dijs[k] <- sqrt(sum((Y[i, ] - Y[j, ]) ^ 2))
      }
      res_dy[i, j] <- stat_fun(dijs)
      res_dy[j, i] <- res_dy[i, j]
    }
  }

  res_dy
}

refine_rep <- function(df, repn, ...) {
  tsmessage("Creating ensemble DY")
  dy_med <- combine_dy(repn, stat_fun = stats::median)
  tsmessage("Creating ensemble Y")
  dy_med_mmds <- smallvis(stats::as.dist(dy_med), method = "mmds", scale = FALSE, Y_init = "pca",
                            epoch = 5, eta = 0.0001, verbose = FALSE)
  varargs <- list(...)
  varargs$Y_init <- dy_med_mmds
  do.call(smallvis, lreplace(X = df, varargs))
}

refine_best_rep <- function(df, repn, ...) {
  best <- best_repn(repn)
  varargs <- list(...)
  varargs$Y_init <- best$Y
  do.call(smallvis, lreplace(X = df, varargs))
}

last <- function(x) {
  x[length(x)]
}

izsplot <- function(res_lap, res_spca, name) {
  graphics::plot(res_lap[[name]]$izs, type = "l", xlab = "iter", ylab = "1/Z",
       main = name)
  graphics::lines(res_spca[[name]]$izs, col = "red")
}

izsplot2 <- function(name, ...) {
  l <- list(...)

  ylim <- c(min(l[[1]][[name]]$izs, l[[2]][[name]]$izs,
                l[[3]][[name]]$izs, l[[4]][[name]]$izs),
            max(l[[1]][[name]]$izs, l[[2]][[name]]$izs,
                l[[3]][[name]]$izs, l[[4]][[name]]$izs))
  graphics::plot(l[[1]][[name]]$izs, type = "l", xlab = "iter", ylab = "1/Z",
                 ylim = ylim,
                 main = paste0(name, " perplexity = ", l[[1]][[name]]$perplexity), lwd = 2)
  graphics::lines(l[[2]][[name]]$izs, col = "red", lwd = 2)
  graphics::lines(l[[3]][[name]]$izs, col = "blue", lwd = 2)
  graphics::lines(l[[4]][[name]]$izs, col = "cyan", lwd = 2)
}

izs_last <- function(name, ...) {
  l <- list(...)
  lasts <- c()
  for (i in 1:length(l)) {
    lasts <- c(lasts, last(l[[i]][[name]]$izs))
  }

  lasts
}

# izs_last_tstr("fashion", res_ssep_lap_u100, res_ssep_spca_u100, res_ssep_rand_u100, res_ssep_randx_u100)
izs_last_tstr <- function(name, ...) {
  paste0(signif(izs_last(name, ...), digits = 4), sep = " | ", collapse = "")
}

izs_last2 <- function(res) {
  lasts <- c()
  for (i in 1:length(res)) {
    lasts <- c(lasts, last(res[[i]]$izs))
  }

  lasts
}

izs_last_tstr2 <- function(res) {
  paste0(signif(izs_last2(res), digits = 4), sep = " | ", collapse = "")
}

izsplot_sr <- function(res_sr1k, res_sr2k, res_sr3k, perp) {
  name <- paste0("perp_", perp)
  graphics::plot(res_sr1k[[name]]$izs, type = "l", xlab = "iter", ylab = "1/Z",
       main = paste0("perplexity = ", perp))
  graphics::lines(res_sr2k[[name]]$izs, col = "red")
  graphics::lines(res_sr3k[[name]]$izs, col = "blue")
}

izsplot_sr_perp <- function(res_sr, title) {
  graphics::plot(res_sr$perp_10$izs, type = "l", xlab = "iter", ylab = "1/Z",
       main = title)
  graphics::lines(res_sr$perp_20$izs, col = "red")
  graphics::lines(res_sr$perp_30$izs, col = "orange")
  graphics::lines(res_sr$perp_40$izs, col = "blue")
  graphics::lines(res_sr$perp_50$izs, col = "purple")

}

izs_perp <- function(perps = seq(10, 50, by = 10), ...) {
  res <- list()
  varargs <- list(...)
  if (!is.null(varargs$ret_extra) && methods::is(varargs$ret_extra, "character")) {
    varargs$ret_extra <- unique(c(varargs$ret_extra, "izs"))
  }
  else {
    varargs$ret_extra <- c("izs")
  }
  for (perp in perps) {
    varargs$perplexity <- perp
    res[[paste0("perp_", perp)]] <- do.call(smallvis, varargs)
  }
  res
}

septsne <- function(perplexity, inp_kernel = "gaussian") {
  lreplace(
    tsne(perplexity = perplexity, inp_kernel = inp_kernel),
    init = function(cost, X, eps = .Machine$double.eps, verbose = FALSE,
                    ret_extra = c()) {
      ret_extra = unique(c(ret_extra, "V"))
      cost <- sne_init(cost, X, perplexity = perplexity, kernel = inp_kernel,
                       symmetrize = "symmetric", normalize = TRUE,
                       verbose = verbose, ret_extra = ret_extra)
      P <- cost$P
      cost$plogp <- colSums(P * log((P + eps)))
      cost$eps <- eps

      cost$V <- cost$V / rowSums(cost$V)
      cost$V <- 0.5 * (cost$V + t(cost$V))

      # Due to normalization scheme Vsum = N
      cost$Vsum <- nrow(cost$V)
      cost$invVsum <- 1 / (cost$Vsum)
      N <- nrow(cost$V)

      # dataset dependent
      cost$invZinc <- 2e-4 / 1000
      cost$invZ <- 0
      cost
    },
    gr = function(cost, Y) {
      cost <- cost_update(cost, Y)

      cost$G <- k2g(Y, 4 * cost$W * cost$invVsum * (cost$V - cost$W * cost$invZ * cost$Vsum))
      cost
    },
    update = function(cost, Y) {
      W <- dist2(Y)
      W <- 1 / (1 + W)
      diag(W) <- 0

      cost$W <- W
      cost$invZ <- cost$invZ + cost$invZinc
      cost
    }
  )
}

perp_explore <- function(X, from = 5, to = base::min(nrow(X) - 1, 300),
                         nperps = to - from + 1, scale = "none",
                         verbose = TRUE) {
  if (methods::is(X, "data.frame")) {
    indexes <- which(vapply(X, is.numeric, logical(1)))
    if (verbose) {
      message("Found ", length(indexes), " numeric columns")
    }
    if (length(indexes) == 0) {
      stop("No numeric columns found")
    }
    X <- X[, indexes]
  }

  X <- scale_input(X, scale, verbose = verbose)

  if (verbose) {
    tsmessage("Calculating perplexities from ", from, " to " , to)
  }

  if (is.null(nperps)) {

  }
  perps <- seq(from = from, to = to, length.out = nperps)

  means <- rep(0, nperps)
  medians <- rep(0, nperps)
  mins <- rep(0, nperps)
  maxs <- rep(0, nperps)

  report_every = max(1, base::round(nperps / 10))

  guesses <- NULL
  for (i in 1:length(perps)) {
    perp <- perps[i]

    pres <- x2aff(X, perplexity = perp, tol = 1e-5, kernel = "gaussian",
          guesses = guesses, verbose = FALSE)
    means[i] <- mean(pres$dint)
    medians[i] <- stats::median(pres$dint)
    mins[i] <- min(pres$dint)
    maxs[i] <- max(pres$dint)
    guesses <- pres$beta

    if (verbose && i %% report_every == 0) {
      tsmessage("Calculated ", i, " / ", nperps)
    }
  }

  names(means) <- perps
  names(medians) <- perps
  names(mins) <- perps
  names(maxs) <- perps

  list(
    means = means,
    medians = medians,
    mins = mins,
    maxs = maxs,
    perps = perps
  )
}








# fuzzy_simplicial_set <- function(X) {
#   for i in range(knn_indices.shape[0]):
#     for j in range(n_neighbors):
#     if knn_indices[i, j] == -1:
#     continue  # We didn't get the full knn for i
#   if knn_indices[i, j] == i:
#     val = 0.0
#     elif knn_dists[i, j] - rhos[i] <= 0.0:
#       val = 1.0
#     else:
#       val = np.exp(-((knn_dists[i, j] - rhos[i]) / (sigmas[i] *
#                                                       bandwidth)))
#
#     rows[i * n_neighbors + j] = i
#     cols[i * n_neighbors + j] = knn_indices[i, j]
#     vals[i * n_neighbors + j] = val
#
#     result = scipy.sparse.coo_matrix((vals, (rows, cols)),
#                                      shape=(X.shape[0], X.shape[0]))
#     result.eliminate_zeros()
#
#     transpose = result.transpose()
#
#     prod_matrix = result.multiply(transpose)
#
#     result = set_op_mix_ratio * (result + transpose - prod_matrix) + \
#     (1.0 - set_op_mix_ratio) * prod_matrix
#
#     result.eliminate_zeros()
#
#     return result
# }

lapnorm <- function(A) {
  # D <- diag(colSums(A))
  Dinvsq <- diag(1 / sqrt(colSums(A)))
  # L <- D - A
  # Dinvsq %*% L %*% Dinvsq
  eye(nrow(A)) - Dinvsq %*% A %*% Dinvsq

  # Dsq <- sqrt(colSums(A))
  # Ln <- -t(A / Dsq) / Dsq
  # diag(Ln) <- 1 + diag(Ln)
  # Ln

  # browser()
}

eye <- function(n) {
  diag(1, nrow = n, ncol = n)
}

laplap <- function(A) {
  D <- colSums(A)
  Lr <- -A / D
  diag(Lr) <- 1 + diag(Lr)

  Dsq <- sqrt(D)
  Ln <- -t(A / Dsq) / Dsq
  diag(Ln) <- 1 + diag(Ln)

  w <- eigen(Lr, symmetric = FALSE)
  v <- eigen(Ln, symmetric = TRUE)

  browser()
}

norm_vec <- function(x) {
  x / sqrt(sum(x * x))
}

pBIC <- function(n = NULL, perp = NULL, kl = NULL) {
  if (is.list(n)) {
    perp <- n$perplexity
    kl <- final_cost(n, digits = Inf)
    n <- nrow(n$Y)
  }
  2 * kl + log(n) * (perp / n)
}

pBIC_cb <- function(res) {
  pBIC_res <- pBIC(n = nrow(res$Y), perp = res$perplexity, kl = final_cost(res, digits = Inf))
  paste0(" pBIC = ", formatC(pBIC_res))
}


zca <- function(xm, scale = FALSE) {
  xm <- scale(xm, scale = scale)
  n <- nrow(xm)
  # This implementation does SVD directly on xm
  # Uses La.svd so that the loadings V are already transposed and we
  # only need to untranspose if ZCA is asked for.
  svdx <- La.svd(xm, nu = 0, nv = ncomp)
  dm <- diag(sqrt(n - 1) / (svdx$d[1:ncomp] + epsilon), nrow = ncomp)
  vm <- svdx$vt[1:ncomp, , drop = FALSE]
  wm <- dm %*% vm
  if (zca) {
    wm <- t(vm) %*% wm
  }
  xm %*% t(wm)
}

smooth_knn_distances_old <- function(X, k = 15, tol = 1e-5,
                                 min_k_dist_scale = 1e-3, verbose = FALSE) {
  if (verbose) {
    tsmessage("Commencing smooth kNN distance calibration for k = ", formatC(k))
  }
  x_is_dist <- methods::is(X, "dist")
  if (x_is_dist) {
    D <- X
    n <- attr(D, "Size")

    D <- as.matrix(D)
  }
  else {
    XX <- rowSums(X * X)
    n <- nrow(X)
  }

  P <- matrix(0, n, n)
  sigma <- rep(1, n)
  rhos <- rep(0, n)
  logU <- log2(k)
  # Will contain index of any point where all neighbors have zero distance
  # Hopefully not too many of these
  badi <- c()
  # Running total of distances, to be converted to a mean if badi is non-empty
  sumD <- 0

  for (i in 1:n) {
    sigma_min <- 0
    sigma_max <- Inf

    if (x_is_dist) {
      Di <- D[i, -i]
    }
    else {
      Di <- (XX[i] + XX - 2 * colSums(tcrossprod(X[i, ], X)))[-i]
      Di[Di < 0] <- 0
      Di <- sqrt(Di)
    }
    sumD <- sumD + sum(Di)
    rho <- min(Di[Di > 0])
    Di[Di < rho] <- rho
    thisP <- exp(-(Di - rho) / sigma[i])
    H <- sum(thisP)

    Hdiff <- H - logU
    tries <- 0

    while (abs(Hdiff) > tol && tries < 50) {
      if (Hdiff > 0) {
        sigma_max <- sigma[i]
        sigma[i] <- 0.5 * (sigma[i] + sigma_min)
      }
      else {
        sigma_min <- sigma[i]
        if (is.infinite(sigma_max)) {
          sigma[i] <- 2 * sigma[i]
        } else {
          sigma[i] <- 0.5 * (sigma[i] + sigma_max)
        }
      }

      thisP <- exp(-(Di - rho) / sigma[i])
      H <- sum(thisP)
      Hdiff <- H - logU
      tries <- tries + 1
    }

    if (rho > 0.0) {
      meanDi <- mean(Di)
      if (sigma[i] < min_k_dist_scale * meanDi) {
        sigma[i] <- min_k_dist_scale * meanDi
        thisP <- exp(-(Di - rho) / sigma[i])
      }
    }
    else {
      badi <- c(badi, i)
    }

    P[i, -i] <- thisP
    rhos[i] <- rho
  }

  # Deal with completely pathological cases:
  # NB this should be impossible if considering all other points in the dataset
  # because it implies that all distances are zero. Might become relevant if
  # we only consider distances of the k-neighborhood
  if (length(badi) > 0) {
    meanD <- sumD / (n * n)
    bad_sigma <- min_k_dist_scale * meanD
    bad_P <- exp(-1 / bad_sigma)
    for (i in 1:length(badi)) {
      sigma[badi[i]] <- bad_sigma
      P[badi[i], -badi[i]] <- bad_P
    }
  }

  if (verbose) {
    summarize(sigma, "sigma summary")
  }
  list(P = P, sigma = sigma, rho = rhos)
}

swiss_ss_n <- function(df, perplexity = 40) {
  ns <- seq(from = 100, to = 1000, by = 25)
  all_izs <- rep(0, length(ns))

  for (i in 1:length(ns)) {
    n <- ns[i]
    all_izs[i] <- last(smallvis(df[1:n, ], perplexity = perplexity, Y_init = "spca", eta = 100, max_iter = 50000, epoch = 100, ret_extra = c("izs"))$izs)
  }
  names(all_izs) <- ns
  all_izs
}

ss_n_sample <- function(df, perplexity = 40) {
  ns <- seq(from = 100, to = 1000, by = 25)
  all_izs <- rep(0, length(ns))

  for (i in 1:length(ns)) {
    n <- ns[i]
    all_izs[i] <- last(smallvis(df[sample(nrow(df), n), ], perplexity = perplexity, Y_init = "spca", eta = 100, max_iter = 50000, epoch = 100, ret_extra = c("izs"))$izs)
  }
  names(all_izs) <- ns
  all_izs
}
