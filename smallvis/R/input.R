# Input Preprocess --------------------------------------------------------

# Scale X according to various strategies
scale_input <- function(X, scale, verbose = FALSE) {
  if (is.null(scale)) {
    scale <- "none"
  }
  if (is.logical(scale)) {
    if (scale) {
      scale <- "scale"
    } else {
      scale <- "none"
    }
  }
  scale <- match.arg(tolower(scale), c("none", "scale", "range", "absmax"))

  X <- as.matrix(X)
  if (anyNA(X)) {
    stop("Input data contains NA: missing data is not allowed")
  }
  switch(
    scale,
    range = {
      tsmessage("Range scaling X")
      X <- as.matrix(X)
      X <- X - min(X)
      span <- max(X)
      if (span == 0) {
        stop("Range scaling requires input data with non-zero range")
      }
      X <- X / span
    },
    absmax = {
      tsmessage("Normalizing by abs-max")
      X <- base::scale(X, scale = FALSE)
      magnitude <- max(abs(X))
      if (magnitude == 0) {
        stop("Abs-max scaling requires input data with non-zero variation")
      }
      X <- X / magnitude
    },
    scale = {
      tsmessage("Scaling to zero mean and unit variance")
      variances <- apply(X, 2, stats::var)
      X <- X[, variances > 0, drop = FALSE]
      if (ncol(X) == 0) {
        stop("Scaling requires at least one non-zero-variance column")
      }
      tsmessage("Kept ", ncol(X), " non-zero-variance columns")
      X <- base::scale(X, scale = TRUE)
    },
    none = {
      X <- as.matrix(X)
    }
  )
  X
}

# Reduce input dimensionality via PCA and also optionally whiten data
pca_preprocess <- function(X, pca, whiten, initial_dims, verbose = FALSE) {
  # We won't do PCA if the rank of the input is less than the requested
  # initial dimensionality
  if (pca) {
    pca <- min(nrow(X) - 1, ncol(X)) >= initial_dims
  }
  if (pca) {
    if (whiten) {
      tsmessage(
        "Reducing initial dimensionality with PCA and ",
        "whitening to ",
        initial_dims,
        " dims"
      )
      X <- pca_whiten(
        X = X,
        ncol = initial_dims,
        verbose = verbose
      )
    } else {
      tsmessage(
        "Reducing initial dimensionality with PCA to ",
        initial_dims,
        " dims"
      )
      X <- pca_scores(
        X = X,
        ncol = initial_dims,
        verbose = verbose
      )
    }
  }
  list(X = X, pca = pca)
}
