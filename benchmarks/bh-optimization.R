# Run in a fresh R process for each installed variant.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("Usage: Rscript bh-optimization.R LIBRARY VARIANT N THREADS STEPS REPS")
}
lib <- normalizePath(args[1])
variant <- match.arg(args[2], c("baseline", "shared"))
n <- as.integer(args[3])
threads <- as.integer(args[4])
steps <- as.integer(args[5])
reps <- as.integer(args[6])
stopifnot(n > 90L, threads >= 0L, steps >= 1L, reps >= 1L)

library(smallvis, lib.loc = lib)
ns <- asNamespace("smallvis")
stopifnot(
  normalizePath(system.file(package = "smallvis")) == file.path(lib, "smallvis")
)
dll <- getLoadedDLLs()[["smallvis"]][["path"]]
symbols <- system2("nm", c("-C", shQuote(dll)), stdout = TRUE)
stopifnot(
  any(grepl("_bh_tsne_negative_gradient_single\\(", symbols)) ==
    (variant == "baseline"),
  any(grepl("bh_sne_gradient_cpp\\(", symbols)) == (variant == "shared")
)

set.seed(20261107)
X <- matrix(rnorm(n * 8), ncol = 8)
Y <- X[, 1:2] / 2
cost <- get("bhtsne", ns)(
  perplexity = 30,
  inp_kernel = "nngaussian",
  nn = "approximate",
  n_threads = threads,
  theta = 0.5
)
cost <- get("cost_init", ns)(cost, X, steps)
opt <- get("opt_create", ns)(list("dbd", eta = n / 12))
step <- get("opt_step", ns)
run <- function() {
  state <- list(cost_fn = cost, opt = opt, Y = Y)
  for (i in seq_len(steps)) {
    state <- step(state$opt, state$cost_fn, state$Y, i)
  }
  state$Y
}
invisible(run())
checksum <- unname(tools::md5sum(dll))
for (i in seq_len(reps)) {
  invisible(gc(verbose = FALSE))
  timing <- system.time(out <- run())
  stopifnot(all(is.finite(out)))
  cat(
    paste(
      variant,
      n,
      threads,
      steps,
      length(cost$P@x),
      timing[["elapsed"]],
      timing[["user.self"]] + timing[["sys.self"]],
      sprintf("%.16g", sum(out^2)),
      checksum,
      sep = ","
    ),
    "\n",
    sep = ""
  )
}
