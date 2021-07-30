#' Identifies outliers in univariate time series using the algorithm lookout.
#'
#'
#' This is the time series implementation of lookout.
#' @param x The input univariate time series.
#' @inheritParams lookout
#' @return A lookout object.
#' @seealso \code{\link{lookout}}
#'
#' @examples
#' set.seed(1)
#' x <- arima.sim(list(order = c(1,1,0), ar = 0.8), n = 200)
#' x[50] <- x[50] + 10
#' plot(x)
#' lo <- lookout_ts(x)
#' lo
#' @export lookout_ts
lookout_ts <- function(x, alpha = 0.05){
  u <- c(0,diff(diff(x)))
  out <- lookout(u, alpha = alpha, unitize = FALSE, bw = NULL, gpd = NULL)

  outliers <- out$outliers[ ,1]

  if(length(outliers) > 1){
    oo <- c()
    clust <- cumsum(c(1, diff(outliers) > 1))
    len <- max(clust)
    for(kk in 1:len){
      inds <- outliers[which(clust==kk)]
      oo <- c(oo, inds[which.min(out$outlier_probability[inds])] )
    }
    inds <- which(out$outliers[ ,1] %in% oo)
    out$outliers <- out$outliers[inds, ]
  }
  out
}
