#' Identifies bandwidth for outlier detection.
#'
#' This function identifies the bandwidth that is used in the kernel density
#' estimate computation. The function uses topological data analysis (TDA)
#' to find the badnwidth.
#'
#' @inheritParams lookout
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
#' find_tda_bw(X, fast = TRUE)
#'
#' @export
find_tda_bw <- function(X, fast = TRUE, gamma = 0.98, use_differences = FALSE) {
  stopifnot(gamma > 0 && gamma <= 1)
  X <- as.matrix(X)

  # select a subset of X for tda computation
  if (fast) {
    inds <- subset_for_tda(X)
    Xsub <- X[inds, ]
  } else {
    Xsub <- X
  }


  # Code above replaced with much faster mlpack computation, which produces identical output
  # Updated with mlpack version fix
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
    m <- NCOL(X)
    return(unname(quantile(death_radi, probs = gamma, type = 8L)))
  }
}
