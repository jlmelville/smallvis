# Piecewise interpolation of the Cauchy kernel and its square, as in
# Linderman et al. (2019), doi:10.1038/s41592-018-0308-4.
fft_tsne <- function(
  perplexity,
  inp_kernel = "nngaussian",
  nn = "approximate",
  eps = .Machine$double.eps,
  n_threads = 0,
  fft_backend = "stats",
  grid_interval = 1
) {
  if (
    !is.character(inp_kernel) ||
      length(inp_kernel) != 1L ||
      !inp_kernel %in% c("nngaussian", "knn")
  ) {
    stop("FFT t-SNE requires inp_kernel = 'nngaussian' or 'knn'")
  }
  nn <- match.arg(nn, c("approximate", "exact"))
  if (
    !is.numeric(n_threads) ||
      length(n_threads) != 1L ||
      !is.finite(n_threads) ||
      n_threads < 0 ||
      n_threads > .Machine$integer.max ||
      n_threads != trunc(n_threads)
  ) {
    stop("n_threads must be a nonnegative integer")
  }
  if (
    !is.numeric(eps) ||
      is.complex(eps) ||
      length(eps) != 1L ||
      !is.finite(eps) ||
      eps < 0
  ) {
    stop("FFT t-SNE requires finite nonnegative eps")
  }
  fft <- fft_settings(fft_backend, grid_interval)
  cost <- sparse_tsne_cost(
    perplexity,
    inp_kernel,
    nn,
    n_threads,
    eps,
    pfn = function(cost, Y) {
      Z <- fft$sums(Y, mass_only = TRUE)$Z
      cost$pcost <- cost$plogp + cauchy_edge_loss(cost$P, Y, log(Z + eps))
      cost
    },
    gr = function(cost, Y) {
      repulsion <- fft$sums(Y)
      P <- cost$P
      attraction <- tsne_sparse_attraction_cpp(P@i, P@p, P@x, Y, n_threads)
      cost$G <- 4 * (attraction - repulsion$force / (repulsion$Z + eps))
      cost
    }
  )
  use_fft_cost(cost, fft)
}

fft_settings <- function(fft_backend, grid_interval) {
  if (
    !is.numeric(grid_interval) ||
      length(grid_interval) != 1L ||
      !is.finite(grid_interval) ||
      grid_interval <= 0
  ) {
    stop("grid_interval must be a finite positive number")
  }
  if (
    !is.character(fft_backend) ||
      length(fft_backend) != 1L ||
      !fft_backend %in% c("stats", "fftwtools")
  ) {
    stop("fft_backend must be 'stats' or 'fftwtools'")
  }
  transform <- stats::fft
  if (fft_backend == "fftwtools") {
    if (!requireNamespace("fftwtools", quietly = TRUE)) {
      stop("Install fftwtools or set fft_backend = 'stats'")
    }
    transform <- function(x, inverse = FALSE) {
      fftwtools::fftw2d(x, inverse = as.integer(inverse))
    }
  }

  list(
    backend = fft_backend,
    sums = function(Y, mass_only = FALSE) {
      fft_tsne_repulsion(Y, grid_interval, transform, mass_only)
    }
  )
}

use_fft_cost <- function(cost, fft) {
  if (!is.null(fft)) {
    cost$repulsion <- "fft"
    cost$fft_backend <- fft$backend
    cost$gradient_only <- TRUE
  }
  cost
}

cauchy_edge_loss <- function(P, Y, offset = 0) {
  col <- rep(seq_len(ncol(P)), diff(P@p))
  row <- P@i + 1L
  dx <- Y[col, 1L] - Y[row, 1L]
  dy <- Y[col, 2L] - Y[row, 2L]
  sum(P@x * (log1p(dx^2 + dy^2) + offset))
}

fft_tsne_repulsion <- function(Y, grid_interval, transform, mass_only = FALSE) {
  xlim <- range(Y[, 1L])
  ylim <- range(Y[, 2L])
  limits <- range(xlim, ylim)
  if (any(!is.finite(limits))) {
    stop("FFT t-SNE coordinates must be finite")
  }
  n <- nrow(Y)
  if (xlim[1L] == xlim[2L] && ylim[1L] == ylim[2L]) {
    return(list(Z = as.double(n) * (n - 1L), force = matrix(0, n, 2L)))
  }

  span <- diff(limits)
  # Match FIt-SNE's square grid and round small grids to fast FFT sizes.
  boxes <- floor(max(50, span / grid_interval))
  fast_boxes <- c(
    50,
    55,
    60,
    65,
    70,
    75,
    80,
    85,
    90,
    96,
    100,
    110,
    120,
    130,
    140,
    150,
    175,
    200
  )
  if (boxes < max(fast_boxes)) {
    boxes <- min(fast_boxes[fast_boxes >= boxes])
  }
  side <- 6 * boxes
  grid_size <- side^2
  if (!is.finite(grid_size) || grid_size > 4e6) {
    stop(
      "FFT t-SNE padded grid (",
      format(side),
      " x ",
      format(side),
      ") exceeds four million cells. Increase grid_interval or use Barnes-Hut."
    )
  }
  dims <- rep(side, 2L)
  width <- span / boxes
  local_basis <- function(values) {
    position <- ((values - limits[1L]) / span) * boxes
    cell <- pmin(floor(position), boxes - 1L)
    # The nodes are at 1/6, 1/2, and 5/6 of a cell.
    u <- 3 * (position - cell) - 1.5
    list(
      start = 3 * cell,
      weight = cbind(u * (u - 1) / 2, 1 - u^2, u * (u + 1) / 2)
    )
  }
  bx <- local_basis(Y[, 1L])
  by <- local_basis(Y[, 2L])
  sx <- rep(1:3, times = 3L)
  sy <- rep(1:3, each = 3L)
  grid_index <- outer(bx$start, sx, "+") + side * (outer(by$start, sy, "+") - 1)
  S <- Matrix::sparseMatrix(
    i = as.vector(t(grid_index)),
    p = seq.int(0L, 9L * n, by = 9L),
    x = as.vector(t(bx$weight[, sx] * by$weight[, sy])),
    dims = c(grid_size, n)
  )

  # The signed offsets on the padded grid implement a linear convolution.
  offset <- seq_len(side) - 1L
  offset[offset >= side / 2] <- offset[offset >= side / 2] - side
  d2 <- (offset * (width / 3))^2
  w <- 1 / (1 + outer(d2, d2, "+"))
  convolve <- function(field, spectrum) {
    spatial <- transform(
      transform(array(field, dims)) * spectrum,
      inverse = TRUE
    )
    as.vector(Re(spatial) / grid_size)
  }

  XY <- Y - (limits[1L] + span / 2)
  charge <- if (mass_only) matrix(1, n, 1L) else cbind(1, XY)
  grid_charge <- as.matrix(S %*% charge)
  mass <- convolve(grid_charge[, 1L], transform(w))
  Z <- sum(Matrix::crossprod(S, mass)) - n
  if (!is.finite(Z) || Z <= 0) {
    stop(
      "FFT t-SNE could not estimate a positive normalization; reduce grid_interval"
    )
  }
  if (mass_only) {
    return(list(Z = Z))
  }

  spectrum <- transform(w^2)
  grid_force <- vapply(
    seq_len(3L),
    function(j) {
      convolve(grid_charge[, j], spectrum)
    },
    numeric(grid_size)
  )
  phi <- as.matrix(Matrix::crossprod(S, grid_force))
  list(Z = Z, force = XY * phi[, 1L] - phi[, 2:3])
}
