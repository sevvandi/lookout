#' Identifies outliers in univariate time series using the algorithm lookout.
#'
#' This is the time series implementation of lookout which identifies outliers
#' in the double differenced time series.
#' @param x The input univariate time series.
#' @inheritParams lookout
#' @param ... Other arguments are passed to \code{\link{lookout}}.
#' @return A lookout object.
#' @seealso \code{\link{lookout}}
#'
#' @examples
#' set.seed(1)
#' x <- arima.sim(list(order = c(1, 1, 0), ar = 0.8), n = 200)
#' x[50] <- x[50] + 10
#' plot(x)
#' lo <- lookout_ts(x)
#' lo
#' @export lookout_ts
lookout_ts <- function(x, scale = FALSE, ...) {
  u <- c(0, diff(diff(x)))
  out <- lookout(u, scale = scale, ...)
  outliers <- out$outliers[, 1]
  # Keep only the most extreme outlier(s) in each consecutive sequence of outliers
  if (length(outliers) > 1) {
    oo <- c()
    clust <- cumsum(c(1, diff(outliers) > 1))
    len <- max(clust)
    for (kk in seq_len(len)) {
      inds <- outliers[which(clust == kk)]
      oo <- c(oo, inds[which.min(out$outlier_probability[inds])])
    }
    inds <- which(out$outliers[, 1] %in% oo)
    out$outliers <- out$outliers[inds, ]
  }
  out
}
