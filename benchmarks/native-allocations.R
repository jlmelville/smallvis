# Developer benchmark; see benchmarks/README.md for setup and measurement limits.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop(
    "Usage: Rscript native-allocations.R LIB VARIANT CASE N THREADS ITERATIONS CSV"
  )
}
lib <- normalizePath(args[[1]], mustWork = TRUE)
variant <- match.arg(args[[2]], c("baseline", "candidate"))
case <- match.arg(args[[3]], c("gradient", "calibration", "embedding"))
n <- as.integer(args[[4]])
threads <- as.integer(args[[5]])
iterations <- as.integer(args[[6]])
stopifnot(n >= 32L, threads >= 1L, iterations >= 1L)
library("smallvis", lib.loc = lib)
ns <- asNamespace("smallvis")
dll <- getLoadedDLLs()[["smallvis"]]
registered <- names(getDLLRegisteredRoutines(dll)$.Call)
removed <- "_smallvis_d2_to_tweight_cpp"
# The old DLL's broken conversion returns only the first n transformed entries.
# Its removal distinguishes the compiled variants without changing the kernels
# being timed. Paths alone would not establish which native code was loaded.
if (variant == "baseline") {
  stopifnot(removed %in% registered)
  probe <- ns$d2_to_tweight_cpp(matrix(1:9, 3), n_threads = 1)
  stopifnot(all(probe[4:9] == 0))
} else {
  stopifnot(!removed %in% registered)
}
set.seed(42)
X <- matrix(stats::rnorm(n * 8), n, 8)
Y <- matrix(stats::rnorm(n * 2, sd = 0.1), n, 2)
if (case == "gradient") {
  # Setup deliberately stays identical between installed variants. Avoid an
  # expensive calibration when measuring only gradient copies and arithmetic.
  W <- matrix(0, n, n)
  for (j in seq_len(n)) {
    W[, j] <- 1 / (1 + (Y[, 1] - Y[j, 1])^2 + (Y[, 2] - Y[j, 2])^2)
  }
  diag(W) <- 0
  P <- matrix(1 / (n * (n - 1)), n, n)
  diag(P) <- 0
  Z <- sum(W)
  run <- function() ns$tsne_grad_cpp(P, W, Z, Y, n_threads = threads)
} else if (case == "calibration") {
  run <- function() ns$find_beta_cpp(X, perplexity = 30, n_threads = threads)
} else {
  run <- function() {
    smallvis(
      X,
      Y_init = Y,
      scale = FALSE,
      perplexity = 30,
      max_iter = iterations,
      epoch = 0,
      use_cpp = TRUE,
      n_threads = threads,
      epoch_callback = FALSE,
      verbose = FALSE
    )
  }
}
# Small warmup avoids warming a full-sized allocation into peak RSS. GNU time's
# process peak still includes package loading and data setup, as documented.
invisible(ns$tsne_grad_cpp(
  matrix(0, 2, 2),
  matrix(1, 2, 2),
  2,
  matrix(0, 2, 2),
  n_threads = threads
))
invisible(gc())
count <- if (case == "embedding") 1L else iterations
start <- proc.time()
for (i in seq_len(count)) {
  result <- run()
  stopifnot(all(is.finite(if (case == "calibration") result$W else result)))
}
elapsed <- (proc.time() - start)[["elapsed"]]
checksum <- sum(if (case == "calibration") result$W else result)
row <- data.frame(
  variant = variant,
  case = case,
  n = n,
  dimensions = 8L,
  threads = threads,
  iterations = iterations,
  seed = 42L,
  elapsed_seconds = elapsed,
  seconds_per_iteration = elapsed / iterations,
  checksum = checksum,
  R = R.version.string,
  platform = R.version$platform,
  library = lib,
  dll = dll[["path"]]
)
utils::write.csv(row, args[[7]], row.names = FALSE)
