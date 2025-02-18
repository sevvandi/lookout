#' Identifies bandwidth for outlier detection.
#'
#' This function identifies the bandwidth that is used in the kernel density
#' estimate computation. The function uses topological data analysis (TDA)
#' to find the badnwidth.
#'
#' @inheritParams lookout
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
find_tda_bw <- function(X, fast = TRUE, bw_para = 0.95, bw_power = 1) {
  stopifnot(bw_para > 0 & bw_para <= 1)
  stopifnot(bw_power > 0 & bw_power <= 1)
  X <- as.matrix(X)

  # select a subset of X for tda computation
  if (fast) {
    inds <- subset_for_tda(X)
    Xsub <- X[inds, ]
  } else {
    Xsub <- X
  }

  if (NCOL(X) == 1L) {
    phom <- TDAstats::calculate_homology(dist(Xsub), format = "distmat")
  } else {
    phom <- TDAstats::calculate_homology(Xsub, dim = 0)
  }
  death_radi <- phom[, 3L]
  # Added so that very small death radi are not chosen
  med_radi <- median(death_radi)
  death_radi_upper <- death_radi[death_radi >= med_radi]
  if (bw_para == 1) {
    dr_thres_diff <- diff(death_radi_upper)
    return(death_radi_upper[which.max(dr_thres_diff)])
  } else {
    prob <- 1 - (1 - bw_para) * length(death_radi_upper)^(bw_power - 1)
    unname(quantile(death_radi_upper, probs = min(1,max(0, prob)), type = 8L))
  }
}
