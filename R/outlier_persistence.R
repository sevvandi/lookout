#' Computes outlier persistence for a range of significance values.
#'
#' This function computes outlier persistence for a range of significance
#' values, using the algorithm lookout, an outlier detection method that uses
#' leave-one-out kernel density estimates and generalized Pareto distributions
#' to find outliers.
#'
#' @param X The input data in a matrix, data.frame, or tibble format. All
#'   columns should be numeric.
#' @param alpha Grid of significance levels.
#' @param st_qq The starting quantile for death radii sequence. This will be
#'   used to compute the starting bandwidth value.
#' @param unitize An option to normalize the data. Default is \code{TRUE},
#'   which normalizes each column to \code{[0,1]}.
#' @param num_steps The length of the bandwidth sequence.
#'
#' @return A list with the following components:
#' \item{\code{out}}{A 3D array of \code{N x num_steps x num_alpha} where
#' \code{N} denotes the number of observations, num_steps denote the length
#'   of the bandwidth sequence and num_alpha denotes the number of significance
#'   levels. This is a binary array and the entries are set to 1 if that
#'   observation is an outlier for that particular bandwidth and significance
#'   level.}
#' \item{\code{bw}}{The set of bandwidth values.}
#' \item{\code{gpdparas}}{The GPD parameters used. }
#' \item{\code{lookoutbw}}{The bandwidth chosen by the algorithm \code{lookout}
#'   using persistent homology.}
#'
#' @examples
#' X <- rbind(
#'   data.frame(x = rnorm(500),
#'              y = rnorm(500)),
#'   data.frame(x = rnorm(5, mean = 10, sd = 0.2),
#'              y = rnorm(5, mean = 10, sd = 0.2))
#' )
#' plot(X, pch = 19)
#' outliers <- persisting_outliers(X, unitize = FALSE)
#' outliers
#' autoplot(outliers)
#' @export

persisting_outliers <- function(X, alpha = seq(0.01, 0.1, by=0.01),
                                st_qq = 0.9, unitize = TRUE, num_steps = 20) {
  # Prepare X matrix
  X <- as.matrix(X)
  if (unitize) {
    X <- unitize(X)
  }

  # Calculate persistent homology
  if (NCOL(X) == 1L) {
    phom <- TDAstats::calculate_homology(dist(X), format = "distmat")
  } else {
    phom <- TDAstats::calculate_homology(X, dim = 0)
  }

  # Find bandwiths
  death_radi <- phom[, 3L]
  qq_st <- quantile(death_radi, probs = st_qq)
  qq_en <- max(death_radi) * sqrt(5)
  bw_vals <- seq(qq_st, qq_en, length.out = num_steps)
  q_thres <- quantile(death_radi, probs = 0.5)
  dr_thres <- death_radi[death_radi >= q_thres]
  dr_thres_diff <- diff(dr_thres)
  max_persist_ind <- which.max(dr_thres_diff)
  ind1 <- min(which(death_radi >= q_thres))
  ind <- max_persist_ind + ind1 - 1L
  bw_fixed <- death_radi[ind] * sqrt(5)

  # Find outliers
  lookoutobj1 <- lookout(X, alpha = 0.05, unitize = FALSE, bw = bw_fixed)
  paras <- lookoutobj1$gpd[1:2]
  output <- array(0, dim = c(dim(X)[1], num_steps, length(alpha)))
  for (i in seq_along(bw_vals)) {
    lookoutobj <- lookout(X, alpha = 0.05, unitize = FALSE,
                               bw = bw_vals[i], gpd = paras)
    for (j in seq_along(alpha)) {
      outinds <- which(lookoutobj$outlier_probability < alpha[j])
      output[outinds, i, j] <- 1
    }
  }

  # Return results
  structure(list(
    out = output,
    bw = bw_vals,
    gpdparas = paras,
    lookoutbw = bw_fixed,
    alpha = alpha,
    call = match.call()
  ), class="persistingoutliers")
}
