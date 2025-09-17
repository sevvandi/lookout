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
