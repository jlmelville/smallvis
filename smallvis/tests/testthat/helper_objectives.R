# Independent final-score oracle, also used for Barnes-Hut at theta = 0.
tsne_reference_cost <- function(P, Y) {
  P <- as.matrix(P)
  W <- 1 / (1 + as.matrix(stats::dist(Y))^2)
  diag(W) <- 0
  Q <- W / sum(W)
  positive <- P > 0
  sum(P[positive] * log(P[positive] / Q[positive]))
}

# Evaluate the documented AB-divergence directly on non-self pairs, without
# production power/log helpers or cached column sums. Fixtures avoid clipping.
absne_reference_cost <- function(P, Y, alpha, lambda) {
  W <- 1 / (1 + as.matrix(stats::dist(Y))^2)
  diag(W) <- 0
  Q <- W / sum(W)
  off_diagonal <- row(P) != col(P)
  p <- P[off_diagonal]
  q <- Q[off_diagonal]
  beta <- lambda - alpha
  if (alpha == 0 && beta == 0) {
    return(sum(log(p / q)^2) / 2)
  }
  if (beta == 0) {
    return(
      sum(p^alpha * log(p^alpha / q^alpha) - p^alpha + q^alpha) /
        alpha^2
    )
  }
  if (lambda == 0) {
    return(sum(log(q^alpha / p^alpha) + p^alpha / q^alpha - 1) / alpha^2)
  }
  if (alpha == 0) {
    return(sum(q^beta * log(q^beta / p^beta) - q^beta + p^beta) / beta^2)
  }
  sum(
    -p^alpha * q^beta + alpha / lambda * p^lambda + beta / lambda * q^lambda
  ) /
    (alpha * beta)
}
