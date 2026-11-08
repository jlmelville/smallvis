#' Best t-SNE Result From Multiple Initializations
#'
#' Run t-SNE multiple times from a random initialization, and return the
#' embedding with the lowest cost.
#'
#' This function ignores any value of \code{Y_init} you set, and uses
#' \code{Y_init = "rand"}.
#'
#' Built-in methods prepare input scaling, PCA, and affinities once. All starts
#' share that input graph, including when using approximate neighbors, and
#' each start resets its coordinates and optimizer state. Custom constructors
#' prepare their input separately for each run. Runs are ranked by the sum of
#' \code{costs} at their returned coordinates, independently of epoch history.
#' Per-run \code{time_secs} excludes the shared input preparation.
#'
#' @param nrep Number of repeats.
#' @param keep_all If \code{TRUE}, then the return value is a list of lists,
#' indexed from 1 .. \code{nrep}, with each entry the result from each
#' \code{\link{smallvis}} run. Otherwise just the result with the lowest error
#' is returned.
#' @param ... Arguments to apply to each \code{\link{smallvis}} run.
#' @return The \code{\link{smallvis}} result with the lowest final cost, or
#' if \code{keep_all} is \code{TRUE} all results as a list, indexed as 1 ..
#' \code{nrep}. If \code{ret_extra} is not \code{FALSE}, then the final costs for all
#' \code{nrep} runs are also included in the return value list as a vector
#' called \code{all_costs}. In this case, if \code{keep_all} is \code{TRUE}, then
#' \code{all_costs} appears as an extra item on all results. Additionally,
#' each result will have an extra entry \code{best_rep}, giving the index of the
#' result with the lowest cost.
#' @examples
#' set.seed(42)
#' best <- smallvis_rep(
#'   nrep = 2, X = iris, perplexity = 30, max_iter = 100,
#'   ret_extra = TRUE, epoch_callback = FALSE, n_threads = 1,
#'   verbose = FALSE
#' )
#' range(best$all_costs)
#'
#' \dontrun{
#' # Return best result out of five random initializations
#' tsne_iris_best <- smallvis_rep(
#'   nrep = 5, X = iris, perplexity = 50, method = "tsne",
#'   ret_extra = TRUE
#' )
#' # How much do the costs vary between runs?
#' range(tsne_iris_best$all_costs)
#' # Display best embedding found
#' plot(tsne_iris_best$Y)
#'
#' # Keep all results
#' # First result is in tsne_iris_rep[[1]], second in tsne_iris_rep[[2]] etc.
#' tsne_iris_rep <- smallvis_rep(
#'   nrep = 5, X = iris, perplexity = 50, method = "tsne",
#'   ret_extra = TRUE, keep_all = TRUE
#' )
#' # Index of result with smallest error is in special list item 'best_rep'
#' best_iris <- tsne_iris_rep[[tsne_iris_rep[[1]]$best_rep]]
#' }
#' @export
smallvis_rep <- function(nrep = 10, keep_all = FALSE, ...) {
  if (nrep < 1) {
    stop("nrep must be 1 or greater")
  }
  varargs <- match_smallvis_args(list(...))
  best_res <- NULL
  best_cost <- Inf
  all_costs <- c()
  # Keep requested return type for final result
  ret_extra <- varargs$ret_extra

  # always return extra so we can find the cost
  if (!should_ret_extra(ret_extra)) {
    varargs$ret_extra <- TRUE
  }

  varargs$Y_init <- "rand"
  ret <- list()
  prepared <- NULL
  if (!is.function(varargs$method)) {
    prepared <- do.call(smallvis_run, c(varargs, list(prepare_only = TRUE)))
  }

  for (i in 1:nrep) {
    # If verbose is not explicitly set to FALSE, it's TRUE by default
    if (nnat(varargs$verbose) || is.null(varargs$verbose)) {
      tsmessage("Starting embedding # ", i, " of ", nrep)
    }
    if (is.null(prepared)) {
      res <- do.call(smallvis, varargs)
    } else {
      res <- do.call(smallvis_run, c(varargs, list(prepared = prepared)))
    }

    if (keep_all) {
      ret[[i]] <- res
    }

    final_cost <- sum(res$costs)
    names(final_cost) <- NULL
    all_costs <- c(all_costs, final_cost)

    if (!keep_all) {
      if (final_cost < best_cost) {
        best_cost <- final_cost
        best_res <- res
      }
    }
  }

  if (keep_all) {
    # if keep_all is TRUE and we asked for extra return info
    # also add the final costs and the index of the best result to each result
    if (should_ret_extra(ret_extra)) {
      best_rep <- which.min(all_costs)
      for (i in 1:nrep) {
        ret[[i]]$all_costs <- all_costs
        ret[[i]]$best_rep <- best_rep
      }
    } else {
      # otherwise just store the Y-coordinates of each result
      for (i in 1:nrep) {
        ret[[i]] <- ret[[i]]$Y
      }
    }
  } else {
    # Only keeping one result
    if (should_ret_extra(ret_extra)) {
      # store info about other results on best result list
      best_res$all_costs <- all_costs
      ret <- best_res
    } else {
      ret <- best_res$Y
    }
  }
  ret
}

# If ret_extra is NULL or FALSE, we aren't returning extra info
# If ret_extra is TRUE or a vector, we are returning extra info
should_ret_extra <- function(ret_extra) {
  !is.null(ret_extra) &&
    ((methods::is(ret_extra, "logical") && ret_extra) ||
      methods::is(ret_extra, "character"))
}

#' Dimensionality Reduction With Perplexity Stepping
#'
#' Carry out dimensionality reduction of a (small) dataset using one of a
#' variety of neighbor embedding methods, using a decreasing value of
#' perplexity to avoid bad local minima.
#'
#' This function uses ideas similar to those in the NeRV (Venna et al., 2010)
#' and JSE (Lee et al., 2013), where to avoid local minima, the initial
#' optimization steps use affinities with larger bandwidths (NeRV) or larger
#' perplexity values (JSE). This implementation uses a series of decreasing
#' perplexity values, as in JSE.
#'
#' For details on the arguments that can be passed to the dimensionality
#' reduction routine, see the help text for \code{\link{smallvis}}.
#'
#' To avoid spending too much extra time in perplexity calibrations, the extra
#' perplexities start at the power of 2 closest to, but not greater than,
#' half the dataset size (in terms of number of objects). Further calibrations
#' are then carried out halving the perplexity each time, until the perplexity
#' specified by the user is reached.
#'
#' The number of iterations spent in the larger perplexity values is specified
#' by the \code{step_iter} parameter. This determines the total number
#' of iterations, e.g. if \code{step_iter = 250} and extra optimizations
#' at a perplexity of 1024, 512, 256, 128 and 64 will be carried out, these will
#' run for 50 iterations each. To keep the number of iterations equivalent to
#' that used by a single run of \code{\link{smallvis}}, the value of
#' \code{step_iter} is subtracted from the value of \code{max_iter} before
#' the optimization at the target perplexity is carried out, e.g. if
#' \code{max_iter = 1000} and \code{step_iter = 250}, the final
#' optimization will run for 750 iterations only.
#'
#' Any value of \code{tol}, early or late exaggeration, and their iteration
#' controls is used only with the final optimization. Final-run defaults are
#' evaluated using its remaining iteration budget. An explicit perplexity in a
#' method list defines the target and is temporarily overridden for preliminary
#' runs. If there are no larger perplexities, or \code{step_iter = 0}, the call
#' is equivalent to a direct \code{smallvis} run.
#'
#' With precomputed \code{neighbors}, the graph must cover the entire schedule,
#' including preliminary perplexities. The neighbor budget is checked before
#' any embedding starts.
#'
#' @param step_iter Number of iterations to carry out the perplexity
#'  stepping. Must be < the value of \code{max_iter} when preliminary runs
#'  are needed. Zero skips stepping.
#' @param ... Arguments to be passed to \code{\link{smallvis}}. See 'Details'
#'  for information on which arguments may be modified or ignored during certain
#'  parts of the embedding.
#' @return The result of the final run of \code{\link{smallvis}} at the target
#'   perplexity.
#' @examples
#' embedding <- smallvis_perpstep(
#'   X = iris, perplexity = 30, Y_init = "spca",
#'   step_iter = 25, max_iter = 100, epoch_callback = FALSE,
#'   n_threads = 1, verbose = FALSE
#' )
#' dim(embedding)
#'
#' \dontrun{
#' # t-SNE on the iris with L-BFGS optimization
#' # The 1000 max_iter is split between 250 iterations at perplexity = 64
#' # and then 750 iterations at perplexity = 40.
#' iris_lbfgs_pstep <- smallvis_perpstep(
#'   step_iter = 250, X = iris, scale = FALSE, verbose = TRUE, Y_init = "spca",
#'   ret_extra = c("DX", "DY"), perplexity = 40, max_iter = 1000, opt = list("l-bfgs")
#' )
#' }
#' @export
#' @references
#' Venna, J., Peltonen, J., Nybo, K., Aidos, H., & Kaski, S. (2010).
#' Information retrieval perspective to nonlinear dimensionality reduction for
#' data visualization.
#' \emph{Journal of Machine Learning Research}, \emph{11}, 451-490.
#'
#' Lee, J. A., Renard, E., Bernard, G., Dupont, P., & Verleysen, M. (2013).
#' Type 1 and 2 mixtures of Kullback-Leibler divergences as cost functions in
#' dimensionality reduction based on similarity preservation.
#' \emph{Neurocomputing}, \emph{112}, 92-108.
smallvis_perpstep <- function(step_iter = 250, ...) {
  final_args <- match_smallvis_args(list(...))
  max_iter <- final_args$max_iter
  if (is.null(max_iter)) {
    max_iter <- 1000
  }
  target_perplexity <- final_args$perplexity
  if (is.list(final_args$method) && !is.null(final_args$method$perplexity)) {
    target_perplexity <- final_args$method$perplexity
  }
  if (is.null(target_perplexity)) {
    target_perplexity <- 30
  }

  X <- final_args$X
  n <- if (methods::is(X, "dist")) attr(X, "Size") else nrow(X)
  validate_perplexity(target_perplexity, n)
  if (length(target_perplexity) != 1) {
    stop("Perplexity stepping requires a scalar target perplexity")
  }
  perps <- scale_perps(n = n, target_perp = target_perplexity)
  if (length(perps) > 0 && step_iter > 0) {
    if (max_iter <= step_iter) {
      stop("max_iter must be > step_iter")
    }
    max_iter_step <- max(1, floor(step_iter / length(perps)))
    preliminary <- final_args
    if (!is.null(final_args$neighbors)) {
      preliminary$neighbor_perplexities <- c(perps, target_perplexity)
    }
    preliminary$ret_extra <- FALSE
    preliminary$max_iter <- max_iter_step
    preliminary$tol <- 0
    preliminary$exaggeration_factor <- 1
    preliminary$late_exaggeration_factor <- 1
    preliminary$epoch <- max_iter_step
    if (is.list(preliminary$method)) {
      preliminary$method$perplexity <- NULL
    }
    for (perp in perps) {
      if (nnat(final_args$verbose)) {
        tsmessage(
          "Optimizing at step perplexity ",
          formatC(perp),
          " for ",
          max_iter_step,
          " iterations"
        )
      }
      preliminary$perplexity <- perp
      Y <- do.call(smallvis_run, preliminary)
      preliminary$neighbor_perplexities <- NULL
      preliminary$Y_init <- Y
    }
    final_args$Y_init <- Y
    final_args$max_iter <- max_iter - step_iter
    max_iter <- final_args$max_iter
  }

  if (nnat(final_args$verbose)) {
    tsmessage(
      "Optimizing at target perplexity ",
      formatC(target_perplexity),
      " for ",
      max_iter,
      " iterations"
    )
  }
  do.call(smallvis, final_args)
}

# Match positional and named arguments against the public API before adding
# private runner controls. Keep omitted arguments omitted so dependent defaults
# such as epoch and stop_lying_iter are evaluated for each actual run.
match_smallvis_args <- function(args) {
  call <- as.call(c(list(quote(smallvis)), args))
  as.list(match.call(definition = smallvis, call = call))[-1]
}

# Utility function for perplexity step
scale_perps <- function(n, target_perp) {
  max_perp <- n / 2
  max_perp <- 2^floor(log(max_perp, 2))
  perp <- max_perp
  perps <- numeric(0)

  if (max_perp > target_perp) {
    perps <- c()
    while (perp > target_perp) {
      perps <- c(perps, perp)
      perp <- perp / 2
    }
  }

  perps
}
