# Output Initialization ---------------------------------------------------

# Initialization of the output coordinates
init_out <- function(Y_init,
                     X,
                     n,
                     ndim,
                     pca_preprocessed,
                     verbose = FALSE) {
  switch(Y_init,
    pca = {
      tsmessage("Initializing from PCA scores")
      pca_init(X, ndim, pca_preprocessed, verbose)
    },
    rand = {
      tsmessage("Initializing from random Gaussian")
      matrix(stats::rnorm(ndim * n, sd = 1), n)
    }
  )
}

pca_init <- function(X, ndim, pca_preprocessed, verbose = FALSE) {
  # If we've already done PCA, we can just take the first two columns
  if (pca_preprocessed) {
    X[, 1:2]
  } else {
    pca_scores(X, ncol = ndim, verbose = verbose)
  }
}

# Laplacian Eigenmap (Belkin & Niyogi, 2002)
# Original formulation solves the generalized eigenvalue problem of the
# unnormalized graph Laplacian: L v = lambda D v, where L = D - A
# and uses the bottom eigenvectors v that result
# (ignoring the constant eigenvector associated with the smallest eigenvalue).
#
# This is equivalent to using the top eigenvectors from the usual
# eigendecomposition of a row-normalized Laplacian P = D^-1 A: P v = lambda' v
# so we don't need to depend on an external package for generalized eigenvalues.
# Note that while the eigenvectors are the same, the eigenvalues are
# different: lambda' = 1 - lambda, but we don't use them with Laplacian
# Eigenmaps anyway.
#
# As we only need to calculate the top ndim + 1 eigenvectors (i.e. normally 3)
# it's incredibly wasteful to calculate all of them. Therefore, if the
# RSpectra library is available, we use that instead, which allows for only the
# top eigenvectors to be extracted. Otherwise, use the slower eigen routine.
# A must be symmetric and positive semi definite, but not necessarily
# normalized in any specific way.
laplacian_eigenmap <- function(A,
                               ndim = 2,
                               use_RSpectra = TRUE) {
  # Equivalent to: D <- diag(colSums(A)); M <- solve(D) %*% A
  # This effectively row-normalizes A: colSums is normally faster than rowSums
  # and because A is symmetric, they're equivalent
  M <- A / colSums(A)
  if (use_RSpectra &&
    requireNamespace("RSpectra", quietly = TRUE, warn.conflicts = FALSE)) {
    tsmessage("Using RSpectra for eigenvectors")
    Re(RSpectra::eigs(M, k = ndim + 1)$vectors[, 2:(ndim + 1)])
  } else {
    tsmessage("Using eigen for eigenvectors")
    eigen(M, symmetric = FALSE)$vectors[, 2:(ndim + 1)]
  }
}

# Use a normalized Laplacian. The UMAP approach, taken from version 0.2.1.
normalized_spectral_init <- function(A,
                                     ndim = 2,
                                     use_RSpectra = TRUE) {
  n <- nrow(A)
  # Normalized Laplacian: clear and close to UMAP code, but very slow in R
  # I <- diag(1, nrow = n, ncol = n)
  # D <- diag(1 / sqrt(colSums(A)))
  # L <- I - D %*% A %*% D

  # A lot faster (order of magnitude when n = 1000)
  Dsq <- sqrt(colSums(A))
  L <- -t(A / Dsq) / Dsq
  diag(L) <- 1 + diag(L)

  if (use_RSpectra &&
    requireNamespace("RSpectra", quietly = TRUE, warn.conflicts = FALSE)) {
    tsmessage("Using RSpectra for eigenvectors")
    k <- ndim + 1
    ncv <- max(2 * k + 1, floor(sqrt(n)))
    opt <- list(
      ncv = ncv,
      maxitr = 5 * n,
      tol = 1e-4
    )
    res <- RSpectra::eigs(L,
      k = k,
      which = "SM",
      opt = opt
    )
    vec_indices <- rev(order(res$values, decreasing = TRUE)[1:ndim])
    res <- Re(res$vectors[, vec_indices])
  } else {
    tsmessage("Using eigen for eigenvectors")
    res <- eigen(L, symmetric = FALSE)
    vec_indices <- order(res$values, decreasing = FALSE)[2:(ndim + 1)]
    res <- Re(res$vectors[, vec_indices])
  }
  res
}

is_spectral_init <- function(init) {
  tolower(init) %in% c("laplacian", "normlaplacian")
}

# Rescale embedding so that the standard deviation is the specified value.
# Default gives initialization like t-SNE, but not random.
shrink_coords <- function(X, sdev = 1e-4) {
  scale(X, scale = apply(X, 2, stats::sd) / sdev)
}

# PCA ---------------------------------------------------------------------


# Calculates a matrix containing the first ncol columns of the PCA scores.
# Returns the score matrix unless ret_extra is TRUE, in which case a list
# is returned also containing the eigenvalues
pca_scores <- function(X,
                       ncol = min(dim(X)),
                       verbose = FALSE,
                       ret_extra = FALSE) {
  if (methods::is(X, "dist")) {
    res_mds <- stats::cmdscale(X,
      x.ret = TRUE,
      eig = TRUE,
      k = ncol
    )

    if (ret_extra || verbose) {
      lambda <- res_mds$eig
      varex <- sum(lambda[1:ncol]) / sum(lambda)
      tsmessage(
        "PCA (using classical MDS): ",
        ncol,
        " components explained ",
        formatC(varex * 100),
        "% variance"
      )
    }
    scores <- res_mds$points
  } else {
    if (ncol < 0.5 * min(dim(X))) {
      res <- irlba::prcomp_irlba(
        X,
        n = ncol,
        retx = TRUE,
        center = TRUE,
        scale = FALSE
      )
      scores <- res$x
      ncol <- ncol(res$rotation)
      varex <- sum(res$sdev[1:ncol]^2) / res$totalvar
      tsmessage(
        "PCA (via irlba): ",
        ncol,
        " components explained ",
        formatC(varex * 100),
        "% variance"
      )
    } else {
      X <- scale(X, center = TRUE, scale = FALSE)
      # do SVD on X directly rather than forming covariance matrix
      s <- svd(X, nu = ncol, nv = 0)
      D <- diag(c(s$d[1:ncol]))
      if (verbose || ret_extra) {
        # calculate eigenvalues of covariance matrix from singular values
        lambda <- (s$d^2) / (nrow(X) - 1)
        varex <- sum(lambda[1:ncol]) / sum(lambda)
        tsmessage(
          "PCA (via SVD): ",
          ncol,
          " components explained ",
          formatC(varex * 100),
          "% variance"
        )
      }
      scores <- s$u %*% D
    }
  }

  if (ret_extra) {
    list(scores = scores, lambda = lambda[1:ncol])
  } else {
    scores
  }
}

# Whiten the data by PCA. This both reduces the dimensionality, but also
# scales the scores by the inverse square root of the equivalent eigenvalue
# so that the variance of each column is 1.
pca_whiten <- function(X,
                       ncol = min(dim(X)),
                       eps = 1e-5,
                       verbose = FALSE) {
  pca <- pca_scores(X,
    ncol = ncol,
    verbose = verbose,
    ret_extra = TRUE
  )
  sweep(pca$scores, 2, sqrt(pca$lambda + eps), "/")
}
