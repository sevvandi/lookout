#' Identifies bandwidth for outlier detection.
#'
#' This function identifies the bandwidth that is used in the kernel density
#' estimate computation. The function uses topological data analysis (TDA)
#' to find the badnwidth.
#'
#' @inheritParams lookout
#' @details The value returned is the raw quantile of the minimum spanning tree
#' edge lengths. \code{\link{lookout}} multiplies it by \code{sqrt(m + 4)},
#' where \code{m = NCOL(X)}, to obtain the support radius of the Epanechnikov
#' kernel, so that the kernel has marginal standard deviation equal to this
#' quantile.
#'
#' @param fast Deprecated and ignored. The bandwidth is always computed
#' using all of \code{X}.
#' @param use_differences If TRUE, the bandwidth is set to the lower point
#' of the maximum Rips death radii differences. If FALSE,
#' the gamma quantile of the Rips death radii is used. Default is FALSE.
#'
#' @return The bandwidth
#'
#' @examples
#' X <- rbind(
#'   data.frame(
#'     x = rnorm(500),
#'     y = rnorm(500)
#'   ),
#'   data.frame(
#'     x = rnorm(5, mean = 10, sd = 0.2),
#'     y = rnorm(5, mean = 10, sd = 0.2)
#'   )
#' )
#' find_tda_bw(X)
#'
#' @export
find_tda_bw <- function(X, fast = NULL, gamma = 0.98, use_differences = FALSE) {
  if (!is.null(fast)) {
    warning(
      "`fast` is deprecated and ignored: the bandwidth is always computed using all of `X`."
    )
  }
  stopifnot(gamma > 0 && gamma <= 1)
  X <- as.matrix(X)

  # The minimum spanning tree is computed on all of X. mlpack::emst is quick
  # enough that sub-sampling is unnecessary, and sub-sampling would change the
  # sample size on which the quantile is based, and so change its scaling in n.
  if (packageVersion("mlpack") < "4.8.0") {
    death_radi <- mlpack::emst(X)$output[, 3]
  } else {
    death_radi <- mlpack::emst(X)[, 3]
  }

  # Added so that very small death radi are not chosen
  if (use_differences) {
    med_radi <- median(death_radi)
    death_radi_upper <- death_radi[death_radi >= med_radi]
    dr_thres_diff <- diff(death_radi_upper)
    return(death_radi_upper[which.max(dr_thres_diff)])
  } else {
    return(unname(quantile(death_radi, probs = gamma, type = 8L)))
  }
}
