# Independent final-score oracle, also used for Barnes-Hut at theta = 0.
tsne_reference_cost <- function(P, Y) {
  P <- as.matrix(P)
  W <- 1 / (1 + as.matrix(stats::dist(Y))^2)
  diag(W) <- 0
  Q <- W / sum(W)
  positive <- P > 0
  sum(P[positive] * log(P[positive] / Q[positive]))
}
