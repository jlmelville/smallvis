# Perplexity Calibration --------------------------------------------------

# Calculates the input affinities from X, such that each normalized row of the
# affinity matrix has the specified perplexity (within the supplied tolerance).
# Returns a list containing the affinities, beta values and intrinsic
# dimensionalities.
# NB set default kernel to "exp" to get results closer to the tsne package.
# This differs from the procedure in the t-SNE paper by exponentially weighting
# the distances, rather than the squared distances.
x2aff <- function(X,
                  perplexity = 15,
                  tol = 1e-5,
                  kernel = "gauss",
                  verbose = FALSE,
                  guesses = NULL) {
  x_is_dist <- methods::is(X, "dist")
  if (x_is_dist) {
    D <- X
    n <- attr(D, "Size")

    D <- as.matrix(D)
    if (kernel == "gauss") {
      D <- D * D
    }
  } else {
    XX <- rowSums(X * X)
    n <- nrow(X)
  }

  nperps <- length(perplexity)
  if (nperps > 1 && nperps != n) {
    stop("Must provide one perplexity per point")
  }

  if (!is.null(guesses) && length(guesses) != n) {
    stop("Initial guess vector must match number of observations in X")
  }

  W <- matrix(0, n, n)
  intd <- rep(0, n)
  if (!is.null(guesses)) {
    beta <- guesses
  } else {
    beta <- rep(1, n)
  }
  if (nperps == 1) {
    logU <- log(perplexity)
  } else {
    perps <- perplexity
  }
  bad_perp <- 0

  for (i in 1:n) {
    if (nperps > 1) {
      perplexity <- perps[i]
      logU <- log(perplexity)
    }
    betamin <- -Inf
    betamax <- Inf

    if (x_is_dist) {
      Di <- D[i, -i]
    } else {
      Di <- (XX[i] + XX - 2 * colSums(tcrossprod(X[i, ], X)))[-i]
      Di[Di < 0] <- 0
      if (kernel == "exp") {
        Di <- sqrt(Di)
      }
    }

    # If we haven't been provided with guesses, then try the initialization used
    # for all points in ELKI according to Schubert & Gertz in "Intrinsic
    # t-Stochastic Neighbor Embedding for Visualization and Outlier Detection: A
    # Remedy Against the Curse of Dimensionality?" Using the last optimized beta
    # seems to be better most of the time based on my testing though, so we'll
    # only use it for the first point.
    if (is.null(guesses) && i == 1) {
      beta[1] <- 0.5 * perplexity / mean(Di)
    }

    sres <- shannon(Di, beta[i])
    H <- sres$H
    Wi <- sres$W
    sumWi <- sres$Z

    Hdiff <- H - logU
    tries <- 0
    while (abs(Hdiff) > tol && tries < 50) {
      if (Hdiff > 0) {
        betamin <- beta[i]
        if (is.infinite(betamax)) {
          beta[i] <- beta[i] * 2
        } else {
          beta[i] <- (beta[i] + betamax) / 2
        }
      } else {
        betamax <- beta[i]
        if (is.infinite(betamin)) {
          beta[i] <- beta[i] / 2
        } else {
          beta[i] <- (beta[i] + betamin) / 2
        }
      }

      sres <- shannon(Di, beta[i])
      H <- sres$H
      Wi <- sres$W
      sumWi <- sres$Z

      Hdiff <- H - logU
      tries <- tries + 1
    }
    if (abs(Hdiff) > tol) {
      # Put a weight of 1/perplexity on the perplexity-nearest neighbors
      bad_perp <- bad_perp + 1
      knn_idx <- order(Di, decreasing = FALSE)[1:max(floor(perplexity), 1)]
      knn_idx[knn_idx >= i] <- knn_idx[knn_idx >= i] + 1
      Wi <- rep(0, n)
      Wi[knn_idx] <- 1 / floor(perplexity)

      intd[i] <- 0
      W[i, ] <- Wi
    } else {
      # if we didn't supply estimates for beta manually, then initialize guess for
      # next point with optimized beta for this point: doesn't save many
      # iterations, but why not?
      if (is.null(guesses) && i < n) {
        beta[i + 1] <- beta[i]
      }
      intd[i] <- intd_x2aff(Di, beta[i], Wi, sumWi, H)
      W[i, -i] <- Wi
    }
  }
  sigma <- sqrt(1 / beta)

  if (bad_perp > 0) {
    tsmessage("Warning: ", bad_perp, " perplexity calibrations failed!")
    warning(bad_perp, " perplexity calibrations failed")
  }

  if (verbose) {
    summarize(sigma, "sigma summary", verbose = verbose)
    summarize(intd, "Dint", verbose = verbose)
  }
  list(W = W, beta = beta, dint = intd)
}

# Calculates affinites based on exponential weighting of D2 with beta
# and returns a list containing:
# W, the affinities; Z, the sum of the affinities; H, the Shannon entropy
# This routine relies specifically on input weights being = exp(-beta * D)
# and calculates the Shannon entropy as log(Z) + beta * sum(W * D) / Z
# where Z is the sum of W.
shannon <- function(D2, beta) {
  W <- exp(-D2 * beta)
  Z <- sum(W)

  if (Z == 0) {
    H <- 0
  } else {
    H <- log(Z) + beta * sum(D2 * W) / Z
  }
  list(W = W, Z = Z, H = H)
}

x2aff_sigma <- function(X,
                        sigma = 1e-3,
                        verbose = FALSE,
                        use_cpp = FALSE,
                        n_threads = 1) {
  x_is_dist <- methods::is(X, "dist")
  if (x_is_dist) {
    D <- X

    D <- as.matrix(D)
    D <- D * D
  } else {
    D <- calc_d2(X, use_cpp = use_cpp, n_threads = n_threads)
  }
  beta <- 1 / (sigma * sigma)
  sres <- shannon(D, beta)
  W <- sres$W
  diag(W) <- 0

  list(W = W, beta = beta)
}

# Create a symmetrized distance matrix based on the k-nearest neighbors
# Non-neighbor distances are set to Inf
knn_dist <- function(X, k, n_threads, verbose) {
  if (methods::is(X, "dist")) {
    # If it's already a distance matrix, find k-smallest distance per column
    # (ignoring self-distances of zero) and set everything larger to Inf
    # (potentially more than k finite distances in the event of ties, not
    # going to worry about that)
    D <- as.matrix(X)
    n <- nrow(D)
    if (k > n - 1) {
      stop("k must be not be > n - 1")
    }
    kdists <- Rfast::colnth(D, rep(k + 1, n))
    for (i in 1:n) {
      Di <- D[, i]
      Di[Di > kdists[i]] <- Inf
      D[, i] <- Di
    }
  } else {
    # Find the k-nearest indexes and distances of X, and set the corresponding
    # distance matrix elements
    n <- nrow(X)
    if (k > n - 1) {
      stop("k must be not be > n - 1")
    }
    tsmessage("Finding ", k + 1, " nearest neighbors")
    knn <- rnndescent::brute_force_knn(X, k = k + 1, n_threads = n_threads)
    knn$idx <- knn$idx[, 2:(k + 1)]
    knn$dist <- knn$dist[, 2:(k + 1)]

    D <- matrix(Inf, nrow = n, ncol = n)
    diag(D) <- 0
    for (i in 1:n) {
      D[i, knn$idx[i, ]] <- knn$dist[i, ]
    }
  }

  # symmetrize
  pmin(D, t(D))
}

# Create the knn graph: D[i, j] = 1 if j is one of i's k-nearest neighbors.
# i is NOT considered a neighbor of itself.
# No symmetrization is carried out.
# Used by knnmmds and knn kernel for SNE
knn_graph <- function(X, k, n_threads, verbose) {
  if (methods::is(X, "dist")) {
    D <- as.matrix(X)
    n <- nrow(D)
    if (k > n - 1) {
      stop("k must be not be > n - 1")
    }
    kdists <- Rfast::colnth(D, rep(k + 1, n))
    for (i in 1:n) {
      Di <- D[, i]
      Di[Di > kdists[i]] <- 0
      D[, i] <- 1
    }
  } else {
    # Find the k-nearest indexes and distances of X, and set the corresponding
    # distance matrix elements
    n <- nrow(X)
    if (k > n - 1) {
      stop("k must be not be > n - 1")
    }

    tsmessage("Finding ", k + 1, " nearest neighbors")
    knn <- rnndescent::brute_force_knn(X, k = k + 1, n_threads = n_threads)
    knn$idx <- knn$idx[, 2:(k + 1)]

    D <- matrix(0, nrow = n, ncol = n)
    for (i in 1:n) {
      D[i, knn$idx[i, ]] <- 1
    }
  }
  D
}



# Multiscale perplexities: P is an average over the results of multiple
# perplexities
# as described by de Bodt et al in
# Perplexity-free t-SNE and twice Student tt-SNE (2018)
msp <- function(X,
                perplexities = NULL,
                tol = 1e-5,
                symmetrize = "symmetric",
                row_normalize = TRUE,
                normalize = TRUE,
                verbose = FALSE,
                guesses = NULL) {
  if (methods::is(X, "dist")) {
    n <- attr(X, "Size")
  } else {
    n <- nrow(X)
  }

  if (is.null(perplexities)) {
    perplexities <- idp_perps(n)
  }
  tsmessage(
    "Calculating multi-scale P with perplexities from ",
    formatC(perplexities[1]),
    " to ",
    formatC(last(perplexities))
  )

  res <- NULL
  for (perplexity in perplexities) {
    tsmessage(
      "Commencing calibration for perplexity = ",
      format_perps(perplexity)
    )
    x2a_res <- x2aff(
      X = X,
      perplexity = perplexity,
      tol = tol,
      kernel = "gauss",
      verbose = verbose,
      guesses = guesses
    )
    P <- x2a_res$W
    Q <- scale_affinities(
      P,
      symmetrize = "symmetric",
      row_normalize = TRUE,
      normalize = TRUE
    )
    if (is.null(res)) {
      res$P <- Q
    } else {
      res$P <- res$P + Q
    }
  }

  if (length(perplexities) > 1) {
    res$P <- res$P / length(perplexities)
  }

  if (is.logical(row_normalize)) {
    tsmessage(
      "Effective perplexity of multiscale P approx = ",
      formatC(stats::median(perpp(res$P)))
    )
  }

  res
}

# Use the Intrinsic Dimensionality Perplexity (IDP)
# Scan through the provided perplexities and use the result which maximizes
# the mean correlation dimension (which is an estimate for the intrinsic
# dimensionality). Stops at the first maxmimum found.
idp <- function(X,
                perplexities = NULL,
                tol = 1e-5,
                verbose = FALSE,
                guesses = NULL) {
  if (methods::is(X, "dist")) {
    n <- attr(X, "Size")
  } else {
    n <- nrow(X)
  }

  if (is.null(perplexities)) {
    perplexities <- idp_perps(n)
  }
  if (verbose) {
    tsmessage(
      "Searching for intrinsic dimensionality with perplexities from ",
      formatC(perplexities[1]),
      " to ",
      formatC(last(perplexities))
    )
  }

  corr_dim_max <- -Inf
  idp <- 0
  idp_res <- NULL
  for (perplexity in perplexities) {
    if (verbose) {
      tsmessage(
        "Commencing calibration for perplexity = ",
        format_perps(perplexity)
      )
    }
    x2a_res <- x2aff(
      X = X,
      perplexity = perplexity,
      tol = tol,
      kernel = "gauss",
      verbose = verbose,
      guesses = guesses
    )
    corr_dim <- mean(x2a_res$dint)
    if (corr_dim <= corr_dim_max) {
      break
    } else {
      corr_dim_max <- corr_dim
      idp <- perplexity
      idp_res <- x2a_res
    }
  }
  if (idp <= 0) {
    stop("Unable to find an IDP: all correlation dimensions were -ve")
  }
  if (verbose) {
    tsmessage(
      "Found IDP at perplexity = ",
      formatC(idp),
      " intrinsic dimensionality = ",
      formatC(corr_dim_max)
    )
  }

  idp_res$idp <- idp
  idp_res
}

# Come up with a set of candidate perplexities for finding the Intrinsic
# Dimensionality Perplexity. Use powers of 2 up to around half the data set
# size, or a perplexity of 128, whichever is smaller. Tries to provide a balance
# of coverage of useful perplexities vs time consumption.
idp_perps <- function(n) {
  max_u <- min(128, max(2, ceiling(n / 2)))
  max_uexp <- floor(log2(max_u))
  min_uexp <- min(2, max_uexp)
  2^(min_uexp:max_uexp)
}

# Is the perplexity argument a string or a list with the first element is a
# string?
perp_method <- function(perplexity) {
  method <- ""
  if (is.character(perplexity) || is.list(perplexity)) {
    if (is.list(perplexity)) {
      method <- perplexity[[1]]
    } else {
      method <- perplexity
    }
  }
  tolower(method)
}

# If the user provided a list like ("idp", c(10, 20, 30)), extract the list
# of numbers as the candidate list of perplexities. Otherwise, return NULL
user_idp_perps <- function(perplexity) {
  perplexities <- NULL
  if (is.list(perplexity) && length(perplexity) == 2) {
    perplexities <- perplexity[[2]]
  }
  perplexities
}

# Some symmetrization options mean "actually, no symmetrization please". This
# function returns the ones that will actually produce a symmetric matrix,
# necessary for symmetric methods (e.g. tsne vs asne).
true_symmetrize_options <- function() {
  c("symmetric", "average", "mutual", "umap", "fuzzy")
}

scale_affinities <- function(P,
                             symmetrize = "symmetric",
                             row_normalize = TRUE,
                             normalize = TRUE) {
  # row normalization before anything else
  if (nnat(row_normalize)) {
    if (symmetrize == "rowsymm") {
      P <- 0.5 * (P + Matrix::t(P))
      symmetrize <- "none"
    }
    P <- P / Matrix::rowSums(P)
  } else if (is.numeric(row_normalize)) {
    P <- row_normalize * P / Matrix::rowSums(P)
  }
  # Symmetrize
  P <- switch(symmetrize,
    none = P,
    symmetric = 0.5 * (P + Matrix::t(P)),
    average = 0.5 * (P + Matrix::t(P)),
    mutual = sqrt(P * Matrix::t(P)),
    umap = fuzzy_set_union(P),
    fuzzy = fuzzy_set_union(P),
    stop("unknown symmetrization: ", symmetrize)
  )
  # Normalize
  if (normalize) {
    P <- P / sum(P)
  }
  P
}

sne_init <- function(cost,
                     X,
                     perplexity,
                     kernel = "gaussian",
                     symmetrize = "symmetric",
                     row_normalize = TRUE,
                     normalize = TRUE,
                     n_threads = 0,
                     use_cpp = FALSE,
                     verbose = FALSE,
                     ret_extra = c()) {
  if (tolower(kernel) == "knn") {
    if (is.character(perplexity) || is.list(perplexity)) {
      stop("Can't use intrinsic dimensionality with knn kernel")
    }
    if (length(perplexity) > 1) {
      stop("Can't use multiple perplexities with knn kernel")
    }
    tsmessage("Using knn kernel with k = ", formatC(perplexity))
    P <- knn_graph(X,
      k = perplexity,
      n_threads = n_threads,
      verbose = verbose
    )
    x2ares <- list(W = P)
  } else if (tolower(kernel) == "skd") {
    P <- smooth_knn_distances(
      X,
      k = perplexity,
      tol = 1e-5,
      n_threads = n_threads,
      verbose = verbose
    )$P
    row_normalize <- FALSE
    x2ares <- list(W = P)
  } else if (perp_method(perplexity) == "idp") {
    perplexities <- NULL
    if (is.list(perplexity) && length(perplexity) == 2) {
      perplexities <- perplexity[[2]]
    }

    x2ares <- idp(X,
      perplexities = perplexities,
      tol = 1e-5,
      verbose = verbose
    )
    P <- x2ares$W
    ret_extra <- unique(c(ret_extra, "idp"))
  } else if (perp_method(perplexity) == "multiscale") {
    perplexities <- NULL
    if (is.list(perplexity) && length(perplexity) == 2) {
      perplexities <- perplexity[[2]]
    }

    mspres <- msp(
      X,
      perplexities = perplexities,
      tol = 1e-5,
      symmetrize = symmetrize,
      row_normalize = row_normalize,
      normalize = normalize,
      verbose = verbose
    )
    cost$P <- mspres$P
    return(cost)
  } else if (tolower(kernel) == "sigma") {
    tsmessage("Using fixed sigma = ", formatC(perplexity))
    x2ares <- x2aff_sigma(
      X,
      sigma = perplexity,
      n_threads = n_threads,
      use_cpp = use_cpp,
      verbose = verbose
    )
    P <- x2ares$W
  } else {
    if (!is.numeric(perplexity)) {
      stop("Unknown perplexity method, '", perplexity[[1]], "'")
    }

    # perpnn options:
    # perpnnks: use exact knn with sparse output
    # perpnnas: use approximate knn with sparse output
    # perpnnkd: use exact knn with dense output
    # perpnnad: use approximate knn with dense output
    if (startsWith(kernel, "perpnn")) {
      n <- nrow(X)
      if (perplexity > n - 1) {
        stop("Perplexity too high for number of points")
      }
      k <- 3 * perplexity
      if (k > n - 1) {
        warning(
          "Perplexity probably too high for number of points,",
          " result may not be meaningful"
        )
        k <- n - 1
      }
      if (kernel %in% c("perpnnks", "perpnnkd")) {
        tsmessage(
          "Finding exact nearest neighbors with k = ",
          k,
          " n_threads = ",
          n_threads
        )
        knn <- rnndescent::brute_force_knn(X, k = k + 1, n_threads = n_threads)
      } else {
        tsmessage(
          "Finding approximate nearest neighbors with k = ",
          k,
          " n_threads = ",
          n_threads
        )
        knn <- rnndescent::rnnd_knn(X, k = k + 1, n_threads = n_threads)
      }
      knn_dist <- knn$dist[, 2:(k + 1)]
      knn_idx <- knn$idx[, 2:(k + 1)]

      ret_sparse <- kernel %in% c("perpnnks", "perpnnas")
      tsmessage(
        "Commencing calibration for perplexity = ",
        format_perps(perplexity),
        " n_threads = ",
        n_threads
      )
      P <- find_beta_knn_cpp(
        knn_dist,
        knn_idx,
        perplexity = perplexity,
        n_threads = n_threads,
        ret_sparse = ret_sparse
      )$P
      if (ret_sparse) {
        P <- Matrix::sparseMatrix(
          i = rep(1:n, each = k),
          j = as.vector(t(knn_idx)),
          x = P,
          dims = c(n, n),
          repr = "C"
        )
      }
    } else {
      if (use_cpp) {
        tsmessage(
          "Commencing calibration for perplexity = ",
          format_perps(perplexity),
          " n_threads = ",
          n_threads
        )
        P <- find_beta_cpp(X, perplexity, tol = 1e-5, n_threads = n_threads)$W
      } else {
        tsmessage(
          "Commencing calibration for perplexity = ",
          format_perps(perplexity)
        )
        x2ares <- x2aff(
          X,
          perplexity,
          tol = 1e-5,
          kernel = kernel,
          verbose = verbose
        )
        P <- x2ares$W
      }
    }
  }

  P <- scale_affinities(
    P,
    symmetrize = symmetrize,
    row_normalize = row_normalize,
    normalize = normalize
  )
  cost$P <- P

  if (!methods::is(P, "sparseMatrix")) {
    if (is.logical(row_normalize)) {
      tsmessage(
        "Effective perplexity of P approx = ",
        formatC(stats::median(perpp(P)))
      )
    }

    for (r in unique(tolower(ret_extra))) {
      switch(r,
        v = {
          cost$V <- x2ares$W
        },
        dint = {
          if (!is.null(x2ares$dint)) {
            cost$dint <- x2ares$dint
          }
        },
        beta = {
          if (!is.null(x2ares$beta)) {
            cost$beta <- x2ares$beta
          }
        },
        adegc = {
          cost$adegc <- 0.5 * rowSums(x2ares$W) + colSums(x2ares$W)
        },
        adegin = {
          cost$adegin <- rowSums(x2ares$W)
        },
        adegout = {
          cost$adegout <- colSums(x2ares$W)
        },
        pdeg = {
          cost$pdeg <- colSums(P)
        },
        idp = {
          if (!is.null(x2ares$idp)) {
            cost$idp <- x2ares$idp
          }
        }
      )
    }
  }
  cost
}

# The intrinsic dimensionality associated with a gaussian affinity vector
# Convenient only from in x2aff, where all these values are available
intd_x2aff <- function(D2, beta, W, Z, H, eps = .Machine$double.eps) {
  P <- W / Z
  -2 * beta * sum(D2 * P * (log(P + eps) + H))
}

shannonpr <- function(P, eps = .Machine$double.eps) {
  P <- P / rowSums(P)
  rowSums(-P * log(P + eps))
}

perpp <- function(P) {
  exp(shannonpr(P))
}
