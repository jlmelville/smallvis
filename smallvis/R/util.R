# Utility Functions -------------------------------------------------------

# Create matrix of squared Euclidean distances
# For low dimension, X %*% t(X) seems to a bit faster than tcrossprod(X)
# Small -ve distances are possible
dist2 <- function(X) {
  D2 <- rowSums(X * X)
  D2 + sweep(X %*% t(X) * -2, 2, t(D2), `+`)
}

calc_d2 <- function(X,
                    use_cpp = FALSE,
                    n_threads = 1) {
  if (use_cpp) {
    dist2_cpp(X, n_threads = n_threads)
  } else {
    safe_dist2(X)
  }
}

calc_d <- function(X,
                   use_cpp = FALSE,
                   n_threads = 1) {
  if (use_cpp) {
    dist_cpp(X, n_threads = n_threads)
  } else {
    sqrt(safe_dist2(X))
  }
}

# Squared Euclidean distances, ensuring no small -ve distances can occur
safe_dist2 <- function(X) {
  D2 <- dist2(X)
  D2[D2 < 0] <- 0
  D2
}

calc_d2tweight <- function(D2,
                           use_cpp = FALSE,
                           n_threads = 1) {
  if (use_cpp) {
    d2_to_tweight_cpp(D2, n_threads = n_threads)
  } else {
    1 / (1 + D2)
  }
}

calc_tweight <- function(X,
                         use_cpp = FALSE,
                         n_threads = 1) {
  if (use_cpp) {
    tweight_cpp(X, n_threads = n_threads)
  } else {
    D2 <- calc_d2(X, use_cpp = use_cpp, n_threads = n_threads)
    1 / (1 + D2)
  }
}

# 2-norm of a vector or matrix
norm2 <- function(X) {
  sqrt(sum(X * X))
}

# Simple time stamp
stime <- function() {
  format(Sys.time(), "%T")
}

# message with a time stamp
tsmessage <- function(...,
                      domain = NULL,
                      appendLF = TRUE,
                      force = FALSE,
                      time_stamp = TRUE) {
  verbose <- get0("verbose", envir = sys.parent())

  if (force || (!is.null(verbose) && verbose)) {
    msg <- ""
    if (time_stamp) {
      msg <- paste0(stime(), " ")
    }
    message(msg, ..., domain = domain, appendLF = appendLF)
    utils::flush.console()
  }
}

# merge lists, where anything non-NULL in l is kept
# e.g.
# all(unlist(update_list(list(a = 1, b = 2), list(a = 10, c = 3))) ==
#     unlist(list(a = 1, b = 2, c = 3)))
lmerge <- function(l, l2) {
  for (name in names(l2)) {
    if (is.null(l[[name]])) {
      l[[name]] <- l2[[name]]
    }
  }
  l
}

# replaces the contents of l with the named arguments
# e.g.
# all(lreplace(c(a = 1, b = 2), a = 10, c = 3) ==
#     c(a = 10, b = 2, c = 3))
lreplace <- function(l, ...) {
  varargs <- list(...)
  for (i in names(varargs)) {
    l[[i]] <- varargs[[i]]
  }
  l
}

# relative tolerance between x and y
reltol <- function(x, y) {
  abs(x - y) / min(abs(x), abs(y))
}

# Check if a value is non-null and true
nnat <- function(x) {
  !is.null(x) && is.logical(x) && x
}

# log vector information
summarize <- function(X, msg = "", verbose = FALSE) {
  summary_X <- summary(X, digits = max(3, getOption("digits") - 3))
  tsmessage(msg, ": ", paste(names(summary_X), ":", summary_X, "|", collapse = ""))
}

# Format perplexity as a string. Could be a scalar or a vector. In the latter
# case, just list the first two values and then ellipses
format_perps <- function(perplexity) {
  if (length(perplexity) > 1) {
    paste0(formatC(perplexity[1]), ", ", formatC(perplexity[2]), "...")
  } else {
    formatC(perplexity)
  }
}

# last item of a vector
last <- function(x) {
  x[length(x)]
}

# remove NULL items from a list
remove_nulls <- function(l) {
  l[!vapply(l, is.null, logical(1))]
}
