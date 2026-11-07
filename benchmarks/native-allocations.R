# Developer benchmark; see benchmarks/README.md for setup and measurement limits.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
  stop(
    paste(
      "Usage: Rscript native-allocations.R",
      "LIB VARIANT REVISION CASE N THREADS ITERATIONS CSV"
    )
  )
}
lib <- normalizePath(args[[1]], mustWork = TRUE)
variant <- match.arg(args[[2]], c("baseline", "candidate"))
revision <- args[[3]]
case <- match.arg(
  args[[4]],
  c(
    "gradient",
    "calibration",
    "embedding",
    "mds-gradient",
    "student-weights",
    "mds-embedding"
  )
)
n <- as.integer(args[[5]])
threads <- as.integer(args[[6]])
iterations <- as.integer(args[[7]])
stopifnot(n >= 32L, threads >= 1L, iterations >= 1L)
library("smallvis", lib.loc = lib)
ns <- asNamespace("smallvis")
dll <- getLoadedDLLs()[["smallvis"]]
# Both revisions produce identical values. Check the compiled helper signatures
# to distinguish their buffer implementations, without a test-only native API.
stopifnot(nzchar(Sys.which("nm")))
symbols <- system2("nm", c("-D", "-C", shQuote(dll[["path"]])), stdout = TRUE)
mds_symbol <- grep(" T mmds_grad\\(", symbols, value = TRUE)
student_symbol <- grep(" T tweight\\(", symbols, value = TRUE)
stopifnot(length(mds_symbol) == 1L, length(student_symbol) == 1L)
borrowed <- grepl(
  "mmds_grad(double const*, double const*,",
  mds_symbol,
  fixed = TRUE
)
direct <- grepl("const&, double*,", student_symbol, fixed = TRUE)
stopifnot(borrowed == (variant == "candidate"))
stopifnot(direct == (variant == "candidate"))
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
} else if (case == "mds-gradient") {
  # Column-wise setup avoids a transient dense conversion dominating peak RSS.
  distances <- function(coordinates) {
    result <- matrix(0, n, n)
    for (j in seq_len(n)) {
      delta <- sweep(coordinates, 2, coordinates[j, ])
      result[, j] <- sqrt(rowSums(delta^2))
    }
    result
  }
  R <- distances(X)
  D <- distances(Y)
  run <- function() ns$mmds_grad_cpp(R, D, Y, 0.01, n_threads = threads)
} else if (case == "student-weights") {
  run <- function() ns$tweight_cpp(Y, n_threads = threads)
} else if (case == "calibration") {
  run <- function() ns$find_beta_cpp(X, perplexity = 30, n_threads = threads)
} else {
  run <- function() {
    smallvis(
      X,
      Y_init = Y,
      method = if (case == "mds-embedding") "mmds" else "tsne",
      scale = FALSE,
      perplexity = 30,
      max_iter = iterations,
      epoch = 0,
      tol = 0,
      opt = list(
        "steepd",
        eta = if (case == "mds-embedding") 0.01 / n else 0.1
      ),
      use_cpp = TRUE,
      n_threads = threads,
      ret_extra = TRUE,
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
embedding <- case %in% c("embedding", "mds-embedding")
count <- if (embedding) 1L else iterations
start <- proc.time()
for (i in seq_len(count)) {
  result <- run()
}
elapsed <- (proc.time() - start)[["elapsed"]]
values <- if (case == "calibration") {
  result$W
} else if (embedding) {
  stopifnot(result$iter == iterations)
  result$Y
} else {
  result
}
# Scan columns after timing without another matrix-sized scratch allocation.
stopifnot(all(vapply(
  seq_len(ncol(values)),
  function(j) {
    all(is.finite(values[, j]))
  },
  logical(1)
)))
checksum <- sum(vapply(
  seq_len(ncol(values)),
  function(j) {
    sum(values[, j]^2)
  },
  numeric(1)
))
row <- data.frame(
  variant = variant,
  revision = revision,
  case = case,
  n = n,
  dimensions = 8L,
  threads = threads,
  iterations = iterations,
  seed = 42L,
  elapsed_seconds = elapsed,
  seconds_per_iteration = elapsed / iterations,
  checksum = checksum,
  final_cost = if (embedding) {
    as.numeric(tail(result$itercosts, 1))
  } else {
    NA_real_
  },
  R = R.version.string,
  platform = R.version$platform,
  library = lib,
  dll = dll[["path"]],
  dll_md5 = unname(tools::md5sum(dll[["path"]])),
  mds_symbol = trimws(mds_symbol),
  student_symbol = trimws(student_symbol),
  blas = unname(extSoftVersion()[["BLAS"]]),
  openblas_threads = Sys.getenv("OPENBLAS_NUM_THREADS"),
  omp_threads = Sys.getenv("OMP_NUM_THREADS")
)
utils::write.csv(row, args[[8]], row.names = FALSE)
