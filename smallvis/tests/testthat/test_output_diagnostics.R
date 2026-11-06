test_that("affinity degrees distinguish incoming, outgoing and mean weights", {
  X <- matrix(sin(seq_len(96)^2), 24, 4)
  for (use_cpp in c(FALSE, TRUE)) {
    res <- smallvis(
      X,
      method = "asne",
      perplexity = 4,
      scale = FALSE,
      use_cpp = use_cpp,
      max_iter = 0,
      ret_extra = c("V", "adegin", "adegout", "adegc"),
      epoch_callback = FALSE,
      verbose = FALSE
    )
    incoming <- colSums(res$V)
    outgoing <- rowSums(res$V)
    expect_gt(max(abs(incoming - outgoing)), 0.1)
    expect_equal(res$adegin, incoming)
    expect_equal(res$adegout, outgoing)
    expect_equal(res$adegc, (incoming + outgoing) / 2)
  }
})

test_that("squared MDS reports ordinary input and output distances", {
  X <- matrix(sin(seq_len(96)^2), 24, 4)
  Y <- cbind(sin(1:24), cos(1:24))
  for (input in list(X, stats::dist(X))) {
    res <- smallvis(
      input,
      method = "smmds",
      Y_init = Y,
      scale = FALSE,
      max_iter = 0,
      ret_extra = c("DX", "DY"),
      epoch_callback = FALSE,
      verbose = FALSE
    )
    expect_equal(res$DX, as.matrix(stats::dist(X)), ignore_attr = TRUE)
    expect_equal(res$DY, as.matrix(stats::dist(Y)), ignore_attr = TRUE)
  }
})

test_that("ASNE-family contributions are each observation's conditional KL", {
  X <- matrix(sin(seq_len(96)^2), 24, 4)
  Y <- cbind(sin(1:24), cos(1:24))
  for (method in list(
    "asne",
    list("basne", beta = 0.5),
    list("btasne", beta = 0.5)
  )) {
    res <- smallvis(
      X,
      method = method,
      perplexity = 4,
      Y_init = Y,
      max_iter = 0,
      ret_extra = c("P", "Q"),
      epoch_callback = FALSE,
      verbose = FALSE
    )
    expected <- vapply(
      1:24,
      function(i) {
        keep <- res$P[i, ] > 0
        sum(res$P[i, keep] * log(res$P[i, keep] / res$Q[i, keep]))
      },
      numeric(1)
    )
    expect_equal(res$costs, expected, tolerance = 1e-10)
    expect_true(all(expected >= 0))
    expect_equal(sum(res$costs), sum(expected), tolerance = 1e-10)
  }
})

test_that("Gaussian output weights exclude self from their normalization shifts", {
  X <- matrix(sin(seq_len(96)^2), 24, 4)
  Y <- cbind(seq_len(24) * 100, 0)
  D2 <- as.matrix(stats::dist(Y))^2
  diag(D2) <- Inf
  for (method in c("asne", "ssne")) {
    expected <- matrix(0, 24, 24)
    if (method == "asne") {
      for (i in 1:24) {
        weights <- exp(-D2[i, ] + min(D2[i, ]))
        expected[i, ] <- weights / sum(weights)
      }
    } else {
      weights <- exp(-D2 + min(D2))
      expected <- weights / sum(weights)
    }
    expected <- pmax(expected, .Machine$double.eps)
    diag(expected) <- 0
    res <- smallvis(
      X,
      method = method,
      perplexity = 4,
      Y_init = Y,
      max_iter = 0,
      ret_extra = c("P", "Q"),
      epoch_callback = FALSE,
      verbose = FALSE
    )
    expect_equal(res$Q, expected, tolerance = 1e-12, ignore_attr = TRUE)
    expect_true(all(is.finite(res$costs)))
    reference <- numeric(24)
    for (i in 1:24) {
      keep <- res$P[i, ] > 0
      reference[i] <- sum(
        res$P[i, keep] * log(res$P[i, keep] / expected[i, keep])
      )
    }
    expect_equal(sum(res$costs), sum(reference), tolerance = 1e-10)
    if (method == "asne") {
      expect_equal(rowSums(res$Q), rep(1, 24), tolerance = 1e-12)
    }
  }
})
