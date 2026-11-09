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
    grid_interval = grid_interval,
    transform = transform,
    sums = function(Y, mass_only = FALSE, alpha = 1, log_weighted = FALSE) {
      fft_sne_sums(Y, grid_interval, transform, mass_only, alpha, log_weighted)
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

fft_sne_grid <- function(Y, grid_interval) {
  xlim <- range(Y[, 1L])
  ylim <- range(Y[, 2L])
  limits <- range(xlim, ylim)
  if (any(!is.finite(limits))) {
    stop("FFT coordinates must be finite")
  }
  n <- nrow(Y)
  if (xlim[1L] == xlim[2L] && ylim[1L] == ylim[2L]) {
    return(NULL)
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
      "FFT padded grid (",
      format(side),
      " x ",
      format(side),
      ") exceeds four million cells. Increase grid_interval or use Barnes-Hut."
    )
  }
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
  weights <- bx$weight[, sx] * by$weight[, sy]
  S <- Matrix::sparseMatrix(
    i = as.vector(t(grid_index)),
    p = seq.int(0L, 9L * n, by = 9L),
    x = as.vector(t(weights)),
    dims = c(grid_size, n)
  )

  # The signed offsets on the padded grid implement a linear convolution.
  offset <- seq_len(side) - 1L
  offset[offset >= side / 2] <- offset[offset >= side / 2] - side
  d2 <- (offset * (width / 3))^2
  list(
    XY = Y - (limits[1L] + span / 2),
    S = S,
    weights = weights,
    side = side,
    width = width,
    d2 = outer(d2, d2, "+")
  )
}

# The interpolation stencil also acts on each point's own charge. Remove
# that contribution before normalizing kernels with a small off-diagonal sum.
fft_self_mass <- function(grid, kernel) {
  sx <- rep(1:3, times = 3L)
  sy <- rep(1:3, each = 3L)
  dx <- outer(sx, sx, "-") %% grid$side
  dy <- outer(sy, sy, "-") %% grid$side
  stencil <- matrix(kernel[1 + dx + grid$side * dy], 9L, 9L)
  rowSums((grid$weights %*% stencil) * grid$weights)
}

fft_grid_convolver <- function(grid, transform, mass_only) {
  grid_size <- grid$side^2
  dims <- rep(grid$side, 2L)
  XY <- grid$XY
  charge <- if (mass_only) matrix(1, nrow(XY), 1L) else cbind(1, XY)
  grid_charge <- as.matrix(grid$S %*% charge)
  function(kernel, fields) {
    spectrum <- transform(kernel)
    spatial <- vapply(
      fields,
      function(j) {
        inverse <- transform(
          transform(array(grid_charge[, j], dims)) * spectrum,
          inverse = TRUE
        )
        as.vector(Re(inverse) / grid_size)
      },
      numeric(grid_size)
    )
    as.matrix(Matrix::crossprod(grid$S, spatial))
  }
}

fft_grid_force <- function(grid, phi) {
  grid$XY * phi[, 1L] - phi[, 2:3]
}

fft_sne_sums <- function(
  Y,
  grid_interval,
  transform,
  mass_only = FALSE,
  alpha = 1,
  log_weighted = FALSE
) {
  n <- nrow(Y)
  grid <- fft_sne_grid(Y, grid_interval)
  if (is.null(grid)) {
    return(list(
      Z = as.double(n) * (n - 1L),
      force = matrix(0, n, 2L),
      log_mass = numeric(n),
      log_force = matrix(0, n, 2L)
    ))
  }

  evaluate <- fft_grid_convolver(grid, transform, mass_only)

  if (alpha == 1) {
    w <- 1 / (1 + grid$d2)
  } else {
    log_w <- if (alpha == 0) -grid$d2 else -log1p(alpha * grid$d2) / alpha
    w <- exp(log_w)
  }
  fields <- if (mass_only) 1L else seq_len(3L)
  # The Gaussian force and mass use the same kernel and spectrum.
  mass <- evaluate(w, if (alpha == 0) fields else 1L)
  Z <- if (alpha == 1) {
    sum(mass[, 1L]) - n
  } else {
    sum(mass[, 1L] - fft_self_mass(grid, w))
  }

  # When the off-diagonal mass is tiny relative to the self contribution,
  # the FFT roundoff is amplified by normalization. Use the exact evaluator.
  if (alpha != 1 && (!is.finite(Z) || Z <= n * sqrt(.Machine$double.eps))) {
    return(NULL)
  }
  if (!is.finite(Z) || Z <= 0) {
    stop(
      "FFT could not estimate a positive normalization. Reduce grid_interval."
    )
  }
  result <- list(Z = Z)
  if (!mass_only) {
    phi <- if (alpha == 0) {
      mass
    } else {
      force_weight <- if (alpha == 1) w^2 else w / (1 + alpha * grid$d2)
      evaluate(force_weight, fields)
    }
    result$force <- fft_grid_force(grid, phi)
  }
  if (log_weighted) {
    log_kernel <- w * (if (alpha == 1) log(w) else log_w)
    log_kernel[w == 0] <- 0
    log_sums <- evaluate(log_kernel, fields)
    result$log_mass <- log_sums[, 1L] - fft_self_mass(grid, log_kernel)
    if (!mass_only) {
      result$log_force <- fft_grid_force(grid, log_sums)
    }
  }
  result
}
