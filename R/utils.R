nsphere <- function(n = 3L, r = 1L, npts = 3145L) {
  x <- matrix(1, nrow = npts, ncol = n)
  for (i in seq(n - 1)) {
    j <- n - i
    theta <- seq(from = -pi, to = pi, length.out = npts)
    theta1 <- sample(theta, size = length(theta))
    x[, seq(j)] <- x[, seq(j)] * cos(theta1)
    x[, (j + 1)] <- x[, (j + 1)] * sin(theta1)
  }
  return(x * r)
}

# Unitize each column of X
unitize <- function(X) {
  for (col in seq_len(NCOL(X))) {
    maxcol <- max(X[, col])
    mincol <- min(X[, col])
    if (maxcol != mincol) {
      X[, col] <- (X[, col] - mincol) / (maxcol - mincol)
    }
  }
  X
}

# Rob's function
# Function to transform multivariate data X to have normal margins
transform_normal <- function(X) {
  X <- as.matrix(X)
  for (i in seq(NCOL(X))) {
    X[, i] <- bestNormalize::yeojohnson(X[, i],
      lower = 0, upper = 1, standardize = FALSE
    )$x.t
  }
  X <- weird::mvscale(X, cov = NULL)
  return(X[, ])
}
