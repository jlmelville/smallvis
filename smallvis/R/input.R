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

  switch(scale,
    range = {
      tsmessage("Range scaling X")
      X <- as.matrix(X)
      X <- X - min(X)
      X <- X / max(X)
    },
    absmax = {
      tsmessage("Normalizing by abs-max")
      X <- base::scale(X, scale = FALSE)
      X <- X / abs(max(X))
    },
    scale = {
      tsmessage("Scaling to zero mean and unit variance")
      X <- Filter(stats::var, X)
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
    pca <- min(nrow(X), ncol(X)) >= initial_dims
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
  X
}
