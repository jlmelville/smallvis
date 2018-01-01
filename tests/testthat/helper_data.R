# convert data frame to matrix using numeric columns
x2m <- function(X) {
  m <- as.matrix(X[, which(vapply(X, is.numeric, logical(1)))])
  attr(m, "dimnames") <- NULL
  m
}

iris10 <- x2m(iris[1:10, ])
iris10_Y <- pca_scores(iris10, ncol = 2)
