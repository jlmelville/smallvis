#' Dimensionality Reduction via Neighbor Embedding
#'
#' Carry out dimensionality reduction of a (small) dataset using one of a
#' variety of neighbor embedding methods.
#'
#' Currently supported embedding methods, which can be used as an argument
#' to the \code{method} parameter are:
#' \enumerate{
#'   \item \code{"tsne"} t-Distributed Stochastic Neighbor Embedding
#'   (van der Maaten and Hinton, 2008).
#'   \item \code{"largevis"} the cost function of the LargeVis algorithm
#'   (Tang et al, 2016).
#' }
#'
#' Note that only the cost function is used from these methods in the context
#' of creating the full distance matrix as part of the gradient calculation.
#' None of the approximations or other speed-ups (e.g. Barnes-Hut or approximate
#' nearest neighbors routines) are used.
#'
#' During the optimization, the vizier package
#' (\url{https://www.github.com/jlmelville/vizier}) is used for visualization.
#' To use a custom callback, provide to the \code{epoch_callback} parameter a
#' function with the following signature:
#'
#' function(Y, iter, cost = NULL)
#'
#' where \code{Y} is the matrix of coordinates, \code{iter} is the current
#' iteration and \code{cost} is the current error value, which will be
#' \code{NULL} the first time this function is called (at iteration 0).
#' The function should have no return value, and presumably will call a plot
#' function. See the "Examples" section for the use of a custom callback.
#' Explicitly set \code{epoch_callback} to \code{NULL} or \code{FALSE} to turn
#' this off.
#'
#' @param X Input coordinates or distance matrix.
#' @param k Number of output dimensions for the embedding.
#' @param scale If \code{TRUE}, scale each column to zero mean and unit
#'   variance. Alternatively, you may specify one of the following strings:
#'   \code{"range"}, which range scales the matrix elements between 0 and 1;
#'   \code{"absmax"}, here the columns are mean centered and then the elements
#'   divided by absolute maximum value; \code{"scale"} does the same as using
#'   \code{TRUE}. To use the input data as-is, use \code{FALSE}, \code{NULL}
#'   or \code{"none"}.
#' @param Y_init How to initialize the output coordinates. One of: \code{"rand"},
#'   which initializes from a Gaussian distribution with mean 0 and standard
#'   deviation 1e-4; \code{"pca"}, which uses the first \code{k} scores of the
#'   PCA: columns are centered, but no scaling beyond that which is applied by
#'   the \code{scale} parameter is carried out; \code{"spca"}, which uses the
#'   PCA scores and then scales each score to a standard deviation of 1e-4; or a
#'   matrix can be used to set the coordinates directly. It must have dimensions
#'   \code{n} by \code{k}, where \code{n} is the number of rows in \code{X}.
#' @param perplexity The target perplexity for parameterizing the input
#'   probabilities.
#' @param inp_kernel The input kernel function. Can be either \code{"gauss"}
#'   (the default), or \code{"exp"}, which uses the unsquared distances.
#'   \code{"exp"} is not the usual literature function, but matches the original
#'   rtsne implementation (and it probably doesn't matter very much).
#' @param max_iter Maximum number of iterations in the optimization.
#' @param pca If \code{TRUE}, apply PCA to reduce the dimensionality of
#'   \code{X} before any perplexity calibration, but after apply any scaling
#'   and filtering. The number of principal components to keep is specified by
#'   \code{initial_dims}. You may alternatively set this value to
#'   \code{"whiten"}, in which case \code{X} is also whitened, i.e. the
#'   principal components are scaled by the inverse of the square root of the
#'   equivalent eigenvalues, so that the variance of component is 1.
#' @param initial_dims If carrying out PCA or whitening, the number of
#'   principal components to keep. Must be no greater than the rank of the input
#'   or no PCA or whitening will be carried out.
#' @param method A neighbor embedding method. See "Details".
#' @param min_cost If the cost falls below this value, the optimization will
#'   stop early.
#' @param epoch_callback Function to call after each epoch. See "Details". By
#'   default the current set of coordinates will be plotted. Set to
#'   \code{FALSE} or \code{NULL} to turn this off.
#' @param epoch After every \code{epoch} number of steps, calculates and
#'   displays the cost value and calls \code{epoch_callback}, if supplied.
#' @param momentum Initial momentum value.
#' @param final_momentum Final momentum value.
#' @param mom_switch_iter Iteration at which the momentum will switch from
#'   \code{momentum} to \code{final_momentum}.
#' @param eta Learning rate value.
#' @param min_gain Minimum gradient descent step size.
#' @param exaggeration_factor Numerical value to multiply input probabilities
#'   by, during the early exaggeration phase. Not used if \code{Y_init} is a
#'   matrix. May also provide the string \code{"ls"}, in which case the
#'   dataset-dependent exaggeration technique suggested by Linderman and
#'   Steinerberger (2017) is used.
#' @param stop_lying_iter Iteration at which early exaggeration is turned
#'   off.
#' @param gamma Weighting term for therepulsive versus attractive forces in the
#'   LargeVis cost function. Used only if \code{method = "largevis"}.
#' @param lveps Epsilon used in the LargeVis gradient to prevent division by
#'   zero. Used only if \code{method = "largevis"}.
#' @param ret_extra If \code{TRUE}, return value is a list containing additional
#'   values associated with the t-SNE procedure; otherwise just the output
#'   coordinates. You may also provide a vector of names of potentially large or
#'   expensive-to-calculate values to return, which will be returned in addition
#'   to those value which are returned when this value is \code{TRUE}. See the
#'   \code{Value} section for details.
#' @param verbose If \code{TRUE}, log progress messages to the console.
#' @return If \code{ret_extra} is \code{FALSE}, the embedded output coordinates
#'   as a matrix. Otherwise, a list with the following items:
#' \itemize{
#' \item{\code{Y}} Matrix containing the embedded output coordinates.
#' \item{\code{N}} Number of objects.
#' \item{\code{origD}} Dimensionality of the input data.
#' \item{\code{scale}} Scaling applied to input data, as specified by the
#'   \code{scale} parameter.
#' \item{\code{Y_init}} Initialization type of the output coordinates, as
#'   specified by the \code{Y_init} parameter, or if a matrix was used, this will
#'   contain the string \code{"matrix"}.
#' \item{\code{iter}} Number of iterations the optimization carried out.
#' \item{\code{time_secs}} Time taken for the embedding, in seconds.
#' \item{\code{perplexity}} Target perplexity of the input probabilities, as
#'   specified by the \code{perplexity} parameter.
#' \item{\code{costs}} Embedding error associated with each observation. This is
#'   the sum of the absolute value of each component of the KL cost that the
#'   observation is associated with, so don't expect these to sum to the
#'   reported KL cost.
#' \item{\code{itercosts}} KL cost at each epoch.
#' \item{\code{stop_lying_iter}} Iteration at which early exaggeration is
#'   stopped, as specified by the \code{stop_lying_iter} parameter.
#' \item{\code{mom_switch_iter}} Iteration at which momentum used in
#'   optimization switches from \code{momentum} to \code{final_momentum}, as
#'   specified by the \code{mom_switch_iter} parameter.
#' \item{\code{momentum}} Momentum used in the initial part of the optimization,
#'   as specified by the \code{momentum} parameter.
#' \item{\code{final_momentum}} Momentum used in the second part of the
#'   optimization, as specified by the \code{final_momentum} parameter.
#' \item{\code{eta}} Learning rate, as specified by the \code{eta} parameter.
#' \item{\code{exaggeration_factor}} Multiplier of the input probabilities
#'   during the exaggeration phase. If the Linderman-Steinerberger exaggeration
#'   scheme is used, this value will have the name \code{"ls"}.
#' \item{\code{pca_dims}} If PCA was carried out to reduce the initial
#'   dimensionality of the input, the number of components retained, as
#'   specified by the \code{initial_dims} parameter.
#' \item{\code{whiten_dims}} If PCA whitening was carried out to reduce the
#'   dimensionality of the input, the number of components retained, as
#'   specified by the \code{initial_dims} parameter.
#' }
#' Additionally, if you set \code{ret_extra} to a vector of names, these will
#' be returned in addition to the values given above. These values are optional
#' and must be explicitly asked for, because they are either expensive to
#' calculate, take up a lot of memory, or both. The available optional values
#' are:
#' \itemize{
#' \item{\code{X}} The input data, after filtering and scaling.
#' \item{\code{P}} The input probabilities.
#' \item{\code{Q}} The output probabilities.
#' \item{\code{DX}} Input distance matrix. The same as \code{X} when the input
#'   data is already a distance matrix.
#' \item{\code{DY}} Output coordinate distance matrix.
#' }
#'
#' @examples
#' \dontrun{
#'
#' # tsne is the default. verbose = TRUE logs progress to console
#' # Also automatically uses github vizier package for plotting coordinates
#' # during optimization
#' tsne_iris <- smallvis(iris, perplexity = 50, verbose = TRUE)
#'
#' # Can use a custom epoch_callback for visualization
#' colors = rainbow(length(unique(iris$Species)))
#' names(colors) = unique(iris$Species)
#' ecb = function(x, y) {
#'   plot(x, t = 'n')
#'   text(x, labels = iris$Species, col = colors[iris$Species])
#' }
#' tsne_iris <- smallvis(iris, epoch_callback = ecb, perplexity = 50, verbose = TRUE)
#'
#' # To turn off visualization entirely:
#' tsne_iris <- smallvis(iris, epoch_callback = FALSE, perplexity = 50, verbose = TRUE)
#'
#' # Try the LargeVis cost function, which also requires a gamma parameter to
#' # be specified:
#' largevis_iris <- smallvis(iris, method = "largevis", gamma = 7,
#'                           epoch_callback = ecb, perplexity = 50, verbose = TRUE)
#'
#' # Use the early exaggeration suggested by Linderman and Steinerberger
#' tsne_iris_ls <- smallvis(iris, epoch_callback = ecb, perplexity = 50,
#'                          exaggeration_factor = "ls")
#'
#' # Make embedding deterministic by initializing with scaled PCA scores
#' tsne_iris_spca <- smallvis(iris, epoch_callback = ecb, perplexity = 50,
#'                            exaggeration_factor = "ls", Y_init = "spca")
#'
#' # Return extra details about the embedding
#' tsne_iris_extra <- smallvis(iris, epoch_callback = ecb, perplexity = 50,
#'                             exaggeration_factor = "ls", Y_init = "spca", ret_extra = TRUE)
#'
#' # Return even more details (which can be slow to calculate or take up a lot of memory)
#' tsne_iris_xextra <- smallvis(iris, epoch_callback = ecb, perplexity = 50,
#'                              exaggeration_factor = "ls", Y_init = "spca",
#'                              ret_extra = c("P", "Q", "X", "DX", "DY"))
#'
#' # Reduce initial dimensionality to 3 via PCA
#' # (But you would normally do this with a much larger dataset)
#' tsne_iris_pca <- smallvis(iris, epoch_callback = ecb, perplexity = 50,
#'                           pca = TRUE, initial_dims = 3)
#'
#' # Or use PCA whitening, so all columns of X have variance = 1
#' tsne_iris_whiten <- smallvis(iris, epoch_callback = ecb, perplexity = 50,
#'                              pca = "whiten", initial_dims = 3)
#' }
#' @references
#' Van der Maaten, L., & Hinton, G. (2008).
#' Visualizing data using t-SNE.
#' \emph{Journal of Machine Learning Research}, \emph{9} (2579-2605).
#' \url{http://www.jmlr.org/papers/v9/vandermaaten08a.html}
#'
#' Tang, J., Liu, J., Zhang, M., & Mei, Q. (2016, April).
#' Visualizing large-scale and high-dimensional data.
#' In \emph{Proceedings of the 25th International Conference on World Wide Web}
#' (pp. 287-297).
#' International World Wide Web Conferences Steering Committee.
#' \url{https://arxiv.org/abs/1602.00370}
#'
#' Linderman, G. C., & Steinerberger, S. (2017).
#' Clustering with t-SNE, provably.
#' \emph{arXiv preprint} \emph{arXiv}:1706.02582.
#' \url{https://arxiv.org/abs/1706.02582}
#' @export
smallvis <- function(X, k = 2, scale = "absmax", Y_init = "rand",
                 perplexity = 30, inp_kernel = "gauss", max_iter = 1000,
                 pca = FALSE, initial_dims = 50,
                 method = "tsne",
                 epoch_callback = TRUE, epoch = base::round(max_iter / 10),
                 min_cost = 0,
                 momentum = 0.5, final_momentum = 0.8, mom_switch_iter = 250,
                 eta = 500, min_gain = 0.01,
                 exaggeration_factor = 4, stop_lying_iter = 100,
                 gamma = 7, lveps = 0.1,
                 ret_extra = FALSE,
                 verbose = TRUE) {

  if (is.logical(epoch_callback)) {
    if (epoch_callback) {
      epoch_callback <- make_smallvis_cb(X)
    }
    else {
      epoch_callback <- NULL
    }
  }
  else if (is.function(epoch_callback)) {
    force(epoch_callback)
  }
  method <- match.arg(tolower(method), c("tsne", "largevis"))

  if (class(pca) == "character" && pca == "whiten") {
    pca <- TRUE
    whiten <- TRUE
  }
  else {
    whiten <- FALSE
  }
  if (pca && initial_dims < k) {
    stop("Initial PCA dimensionality must be larger than desired output ",
         "dimension")
  }

  start_time <- NULL
  ret_optionals <- c()
  if (methods::is(ret_extra, "character")) {
    ret_optionals <- ret_extra
    ret_extra <- TRUE
  }

  if (ret_extra) {
    start_time <- Sys.time()
  }

  if (methods::is(X, "dist")) {
    n <- attr(X, "Size")
  }
  else {
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
    X <- pca_preprocess(X, pca, whiten, initial_dims, verbose = verbose)
    n <- nrow(X)
  }

  if (!is.null(Y_init)) {
    if (methods::is(Y_init, "matrix")) {
      if (nrow(Y_init) != n || ncol(Y_init) != k) {
        stop("Y_init matrix does not match necessary configuration for X")
      }
      Y <- Y_init
      Y_init <- "matrix"
      exaggeration_factor <- 1
    }
    else {
      Y_init <- match.arg(tolower(Y_init), c("rand", "pca", "spca"))
      Y <- init_out(Y_init, X, k, pca_preprocessed = pca,
                    verbose = verbose)
    }
  }
  # Display initialization
  if (!is.null(epoch_callback)) {
    do_callback(epoch_callback, Y, 0)
  }

  itercosts <- c()
  if (tolower(exaggeration_factor) == "ls") {
    # Linderman-Steinerberger exaggeration
    exaggeration_factor <- 0.1 * n
    names(exaggeration_factor) <- "ls"
  }
  else {
    names(exaggeration_factor) <- "ex"
  }

  if (max_iter < 1) {
    return(ret_value(Y, ret_extra, X, scale, Y_init, iter = 0,
                     start_time = start_time, optionals = ret_optionals,
                     pca = ifelse(pca && !whiten, initial_dims, 0),
                     whiten = ifelse(pca && whiten, initial_dims, 0)))
  }

  P <- x2p(X, perplexity, tol = 1e-5, kernel = inp_kernel, verbose = verbose)$P
  P <- 0.5 * (P + t(P))
  P <- P / sum(P)

  if (names(exaggeration_factor) == "ls") {
    if (verbose) {
      message("Linderman-Steinerberger exaggeration = ", formatC(exaggeration_factor))
    }
  }
  P <- P * exaggeration_factor

  uY <- matrix(0, n, k)
  gains <- matrix(1, n, k)
  mu <- momentum
  eps <- .Machine$double.eps
  Z <- 0
  for (iter in 1:max_iter) {
    # D2
    W <- dist2(Y)
    # W
    W <- 1 / (1 + W)
    diag(W) <- 0
    # Force constant (aka stiffness)
    if (method == "tsne") {
      Z <- sum(W)
      G <- 4 * W * (P - W / Z)
    }
    else {
      G <- 4 * P * W - ((gamma * W * W) / ((1 - W) + lveps))
    }

    G <- Y * rowSums(G) - (G %*% Y)

    if (names(exaggeration_factor) == "ls" && iter <= stop_lying_iter) {
      # during LS exaggeration, use gradient descent only with eta = 1
      uY <- -G
    }
    else {
      # compare signs of G with -update (== previous G, if we ignore momentum)
      # abs converts TRUE/FALSE to 1/0
      dbd <- abs(sign(G) != sign(uY))
      gains <- (gains + 0.2) * dbd + (gains * 0.8) * (1 - dbd)
      gains[gains < min_gain] <- min_gain
      uY <- mu * uY - eta * gains * G
    }

    # Update
    Y <- Y + uY

    if (iter == mom_switch_iter) {
      mu <- final_momentum
      if (verbose) {
        message("Switching to final momentum ", formatC(final_momentum),
                " at iter ", iter)
      }
    }

    if (iter == stop_lying_iter && Y_init != "matrix") {
      if (verbose) {
        message("Switching off exaggeration at iter ", iter)
      }
      P <- P / exaggeration_factor
    }

    if (iter %% epoch == 0 || iter == max_iter) {
      # Recenter Y during epoch only
      Y <- sweep(Y, 2, colMeans(Y))

      if (method == "tsne") {
        cost <- sum(P * log((P + eps) / ((W / Z) + eps)))
      }
      else {
        cost <- sum(-P * log(W + eps) - gamma * log1p(-W + eps))
      }

      if (verbose) {
        message(stime(), " Iteration #", iter, " error: ",
                formatC(cost)
                , " ||G||2 = ", formatC(sqrt(sum(G * G))))
      }

      if (!is.null(epoch_callback)) {
        do_callback(epoch_callback, Y, iter, cost)
      }

      if (ret_extra) {
        names(cost) <- iter
        itercosts <- c(itercosts, cost)
      }

      if (cost < min_cost) {
        break
      }
    }
  }

  ret_value(Y, ret_extra, X, scale, Y_init, iter, start_time,
            P, W / Z, eps, perplexity, itercosts,
            stop_lying_iter, mom_switch_iter, momentum, final_momentum, eta,
            exaggeration_factor, optionals = ret_optionals,
            pca = ifelse(pca && !whiten, initial_dims, 0),
            whiten = ifelse(pca && whiten, initial_dims, 0))
}

#' Best t-SNE Result From Multiple Initializations
#'
#' Run t-SNE multiple times from a random initialization, and return the
#' embedding with the lowest cost.
#'
#' This function ignores any value of \code{Y_init} you set, and uses
#' \code{Y_init = "rand"}.
#'
#' @param nrep Number of repeats.
#' @param ... Arguments to apply to each \code{\link{smallvis}} run.
#' @return The \code{\link{smallvis}} result with the lowest final cost.
#' If \code{ret_extra} is not \code{FALSE}, then the final costs for all
#' \code{nrep} runs are also included in the return value list as a vector
#' called \code{all_costs}.
#' @examples
#' \dontrun{
#' # Return best result out of five random initializations
#' tsne_iris_best <- smallvis_rep(nrep = 5, iris, perplexity = 50, method = "tsne", ret_extra = TRUE)
#' # How much do the costs vary between runs?
#' range(tsne_iris_best$all_costs)
#' # Display best embedding found
#' plot(tsne_iris_best$Y)
#' }
#' @export
smallvis_rep <- function(nrep = 10, ...) {
  if (nrep < 1) {
    stop("nrep must be 1 or greater")
  }
  varargs <- list(...)
  best_res <- NULL
  best_cost <- Inf
  all_costs <- c()
  # Keep requested return type for final result
  ret_extra <- varargs$ret_extra
  # Inside loop, always return extra, so we can find the cost
  if (!methods::is(ret_extra, "character") && !ret_extra) {
    varargs$ret_extra <- TRUE
  }

  varargs$Y_init <- "rand"
  for (i in 1:nrep) {
    if (!is.null(varargs$verbose) && varargs$verbose) {
      message(stime(), " Starting embedding # ", i, " of ", nrep)
    }
    res <- do.call(smallvis, varargs)
    final_cost <- res$itercosts[length(res$itercosts)]

    if (final_cost < best_cost) {
      best_cost <- final_cost
      best_res <- res
    }
    names(final_cost) <- NULL
    all_costs <- c(all_costs, final_cost)
  }

  if (!methods::is(ret_extra, "character") && !ret_extra) {
    best_res <- best_res$Y
  }
  else {
    best_res$all_costs <- all_costs
  }
  best_res
}

# Input Preprocess --------------------------------------------------------

# Scale X according to various strategies
scale_input <- function(X, scale, verbose = FALSE) {
  if (is.null(scale)) {
    scale <- "none"
  }
  if (is.logical(scale)) {
    if (scale) {
      scale <- "scale"
    }
    else {
      scale <- "none"
    }
  }
  scale <- match.arg(tolower(scale), c("none", "scale", "range", "absmax"))

  switch(scale,
         range = {
           if (verbose) {
             message(stime(), " Range scaling X")
           }
           X <- as.matrix(X)
           X <- X - min(X)
           X <- X / max(X)
         },
         absmax = {
           if (verbose) {
             message(stime(), " Normalizing by abs-max")
           }
           X <- base::scale(X, scale = FALSE)
           X <- X / abs(max(X))
         },
         scale = {
           if (verbose) {
             message(stime(), " Scaling to zero mean and unit variance")
           }
           X <- Filter(stats::var, X)
           if (verbose) {
             message("Kept ", ncol(X), " non-zero-variance columns")
           }
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
      if (verbose) {
        message(stime(), " Reducing initial dimensionality with PCA and ",
                "whitening to ", initial_dims, " dims")
      }
      X <- pca_whiten(X = X, ncol = initial_dims, verbose = verbose)
    }
    else {
      if (verbose) {
        message(stime(), " Reducing initial dimensionality with PCA to ",
                initial_dims, " dims")
      }
      X <- pca_scores(X = X, ncol = initial_dims, verbose = verbose)
    }
  }
  X
}

# Output Initialization ---------------------------------------------------

# Initialization of the output coordinates
init_out <- function(Y_init, X, ndim, pca_preprocessed, verbose = FALSE) {
  n <- nrow(X)
  switch(Y_init,
         pca = {
           if (verbose) {
             message(stime(), " Initializing from PCA scores")
           }
           if (pca_preprocessed) {
             X[, 1:2]
           }
           else {
             pca_scores(X, ncol = ndim, verbose = verbose)
           }
         },
         spca = {
           if (verbose) {
             message(stime(), " Initializing from scaled PCA scores")
           }
           # If we've already done PCA, we can just take the first two columns
           if (pca_preprocessed) {
             scores <- X[, 1:2]
           }
           else {
             scores <- pca_scores(X, ncol = ndim, verbose = verbose)
           }
           scale(scores, scale = apply(scores, 2, stats::sd) / 1e-4)
         },
         rand = {
           if (verbose) {
             message(stime(), " Initializing from random Gaussian with sd = 1e-4")
           }
           matrix(stats::rnorm(ndim * n, sd = 1e-4), n)
         }
  )
}

# Epoch Functions ---------------------------------------------------------

# Helper function for epoch callback, allowing user to supply callbacks with
# multiple arities.
do_callback <- function(cb, Y, iter, cost = NULL) {
  nfs <- length(formals(cb))
  if (nfs == 1) {
    cb(Y)
  }
  else if (nfs == 2) {
    cb(Y, iter)
  }
  else if (nfs == 3) {
    cb(Y, iter, cost)
  }
}

# Create a callback for visualization
make_smallvis_cb <- function(df) {
  force(df)
  function(Y, iter, cost = NULL) {
    title <- paste0("iter: ", iter)
    if (!is.null(cost)) {
      title <- paste0(title, " cost = ", formatC(cost))
    }
    vizier::embed_plot(Y, df, title = title)
  }
}

# Result Export -----------------------------------------------------------

# Prepare the return value.
# If ret_extra is TRUE, return a list with lots of extra info.
# Otherwise, Y is returned directly.
# If ret_extra is TRUE and iter > 0, then all the NULL-default parameters are
# expected to be present. If iter == 0 then the return list will contain only
# scaling and initialization information.
ret_value <- function(Y, ret_extra, X, scale, Y_init, iter, start_time = NULL,
                      P = NULL, Q = NULL,
                      eps = NULL, perplexity = NULL, pca = 0, whiten = 0,
                      itercosts = NULL,
                      stop_lying_iter = NULL, mom_switch_iter = NULL,
                      momentum = NULL, final_momentum = NULL, eta = NULL,
                      exaggeration_factor = NULL, optionals = c()) {
  if (ret_extra) {
    end_time <- Sys.time()

    if (methods::is(X, "dist")) {
      N <- attr(X, "Size")
      origD <- NULL
    }
    else {
      N <- nrow(X)
      origD <- ncol(X)
    }

    res <- list(
      Y = Y,
      N = N,
      origD = origD,
      scale = scale,
      Y_init = Y_init,
      iter = iter,
      time_secs = as.numeric(end_time - start_time, units = "secs")
    )

    if (pca > 0) {
      res$pca_dims <- pca
    }
    else if (whiten > 0) {
      res$whiten_dims <- whiten
    }

    if (iter > 0) {
      # this was already calculated in the final epoch, but unlikely to be
      # worth the extra coupling and complication of getting it over here
      costs <- colSums(P * log((P + eps) / (Q + eps)))

      if (names(exaggeration_factor) != "ls") {
        names(exaggeration_factor) <- NULL
      }

      res <- c(res, list(
        perplexity = perplexity,
        costs = costs,
        itercosts = itercosts,
        stop_lying_iter = stop_lying_iter,
        mom_switch_iter = mom_switch_iter,
        momentum = momentum,
        final_momentum = final_momentum,
        eta = eta,
        exaggeration_factor = exaggeration_factor
      ))
    }

    for (o in tolower(unique(optionals))) {
      if (o == "p") {
        if (!is.null(P)) {
          res$P <- P
        }
      }
      else if (o == "q") {
        if (!is.null(Q)) {
          res$Q <- Q
        }
      }
      else if (o == "x") {
        res$X <- X
      }
      else if (o == "dx") {
        if (methods::is(X, "dist")) {
          res$DX <- X
        }
        else {
          res$DX <- dist2(X)
        }
      }
      else if (o == "dy") {
        res$DY <- dist2(Y)
      }
    }

    res
  }
  else {
    Y
  }
}

# PCA ---------------------------------------------------------------------


# Calculates a matrix containing the first ncol columns of the PCA scores.
# Returns the score matrix unless ret_extra is TRUE, in which case a list
# is returned also containing the eigenvalues
pca_scores <- function(X, ncol = min(dim(X)), verbose = FALSE,
                       ret_extra = FALSE) {
  X <- scale(X, center = TRUE, scale = FALSE)
  # do SVD on X directly rather than forming covariance matrix
  ncomp <- ncol
  s <- svd(X, nu = ncomp, nv = 0)
  D <- diag(c(s$d[1:ncomp]))
  if (verbose || ret_extra) {
    # calculate eigenvalues of covariance matrix from singular values
    lambda <- (s$d ^ 2) / (nrow(X) - 1)
    varex <- sum(lambda[1:ncomp]) / sum(lambda)
    message("PCA: ", ncomp, " components explained ", formatC(varex * 100),
            "% variance")
  }
  scores <- s$u %*% D
  if (ret_extra) {
    list(
      scores = scores,
      lambda = lambda[1:ncomp]
    )
  }
  else {
    scores
  }
}

# Whiten the data by PCA. This both reduces the dimensionality, but also
# scales the scores by the inverse square root of the equivalent eigenvalue
# so that the variance of each column is 1.
pca_whiten <- function(X, ncol = min(dim(X)), eps = 1e-5, verbose = FALSE) {
  pca <- pca_scores(X, ncol = ncol, verbose = verbose, ret_extra = TRUE)
  sweep(pca$scores, 2, sqrt(pca$lambda + eps), "/")
}


# Perplexity Calibration --------------------------------------------------

# Calculates the input probabilities from X, such that each row probability
# distribution has the specified perplexity (within the supplied tolerance).
# Returns a list containing the probabilities and beta values.
# NB the default kernel, "exp", differs from the procedure in the TSNE paper by
# exponentially weighting the distances, rather than the squared distances.
# Set the kernel to "gauss" to get the squared distance version.
x2p <- function(X, perplexity = 15, tol = 1e-5, kernel = "exp",
                verbose = FALSE) {
  x_is_dist <- methods::is(X, "dist")
  if (x_is_dist) {
    D <- X
    n <- attr(D, "Size")

    D <- as.matrix(D)
    if (kernel == "gauss") {
      D <- D * D
    }
  }
  else {
    XX <- rowSums(X * X)
    n <- nrow(X)
  }

  P <- matrix(0, n, n)
  beta <- rep(1, n)
  logU <- log(perplexity)

  for (i in 1:n) {
    betamin <- -Inf
    betamax <- Inf

    if (x_is_dist) {
      Di <- D[i, -i]
    }
    else {
      Di <- (XX[i] + XX - 2 * colSums(tcrossprod(X[i, ], X)))[-i]
      Di[Di < 0] <- 0
      if (kernel == "exp") {
        Di <- sqrt(Di)
      }
    }
    # Initialization used for all points in ELKI according to Schubert & Gertz
    # in "Intrinsic t-Stochastic Neighbor Embedding for Visualization and
    # Outlier Detection: A Remedy Against the Curse of Dimensionality?"
    # Using the last optimized beta seems to be better most of the time based
    # on my testing though, so we'll only use it for the first point.
    if (i == 1) {
      beta[1] <- 0.5 * perplexity / mean(Di)
    }

    hbeta <- dist_to_prob(Di, beta[i])
    H <- hbeta$H
    thisP <- hbeta$P

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

      hbeta <- dist_to_prob(Di, beta[i])
      H <- hbeta$H
      thisP <- hbeta$P
      Hdiff <- H - logU
      tries <- tries + 1
    }
    # initialize guess for next point with optimized beta for this point
    # doesn't save many iterations, but why not?
    if (i < n) {
      beta[i + 1] <- beta[i]
    }
    P[i, -i] <- thisP
  }
  sigma <- sqrt(1 / beta)

  if (verbose) {
    summary_sigma <- summary(sigma, digits = max(3, getOption("digits") - 3))
    message(stime(), " sigma summary: ",
            paste(names(summary_sigma), ":", summary_sigma, "|", collapse = ""))
  }
  list(P = P, beta = beta)
}

# Given a vector of squared distances and an exponential parameter beta,
# calculates the probabilities and corresponding Shannon entropy.
#
# This routine relies specifically on input weights being = exp(-beta * D)
# and calculates the Shannon entropy as log(Z) + beta * sum(W * D) / Z
# where Z is the sum of W.
#
# Returns a list containing the Shannon entropy and the probability.
dist_to_prob <- function(D, beta) {
  P <- exp(-D * beta)
  Z <- sum(P)
  if (Z == 0) {
    H <- 0
    P <- D * 0
  }
  else {
    H <- log(Z) + beta * sum(D * P) / Z
    P <- P / Z
  }
  list(H = H, P = P)
}

# Utility Functions -------------------------------------------------------

# Create matrix of squared Euclidean distances
# For low dimension, X %*% t(X) seems to a bit faster than tcrossprod(X)
dist2 <- function(X) {
  D2 <- rowSums(X * X)
  D2 <- D2 + sweep(X %*% t(X) * -2, 2, t(D2), `+`)
  D2[D2 < 0] <- 0
  D2
}

# 2-norm of a vector or matrix
norm2 <- function(X) {
  sqrt(sum(X * X))
}

# Simple time stamp
stime <- function() {
  format(Sys.time(), "%T")
}
