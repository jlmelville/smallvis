# convert data frame to matrix using numeric columns
x2m <- function(X) {
  if (!methods::is(X, "matrix")) {
    m <- as.matrix(X[, which(vapply(X, is.numeric, logical(1)))])
    attr(m, "dimnames") <- NULL
  }
  else {
    m <- X
  }
  m
}

# convert dataframe to distance matrix
x2d <- function(X) {
  sqrt(safe_dist2(x2m(X)))
}

iris10 <- x2m(iris[1:10, ])
iris10_Y <- pca_scores(iris10, ncol = 2)

# 10 iris entries which are unique, otherwise knn tests are hard to get
# unique indices and ordering
uiris <- unique(iris)
uirism <- as.matrix(uiris[, -5])
ui10 <- uirism[6:15, ]