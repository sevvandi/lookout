#' Computes outlier persistence for a range of significance values.
#'
#' This function computes outlier persistence for a range of significance values, using the algorithm lookout, an outlier detection method that uses leave-one-out kernel density estimates and generalized Pareto distributions to find outliers.
#'
#' @param X The input data in a dataframe, matrix or tibble format.
#' @param alpha_min The smallest level of significance, set to \code{0.01}.
#' @param alpha_max The maximum level of significance, set to \code{0.1}.
#' @param step The step size, set to \code{0.01}.
#' @param st_qq The starting quantile for death radii sequence. This will be used to compute the starting bandwidth value.
#' @param en_sc The scaling for the maximum death radii. This will be used to compute the biggest bandwidth value.
#' @param unitize An option to normalize the data. Default is \code{TRUE}, which normalizes each column to \code{[0,1]}.
#' @param num_steps The length of the bandwidth sequence.
#'
#' @return A list with the following components:
#' \item{\code{out}}{A 3D array of \code{N x num_steps x num_alpha} where \code{N} denotes the number of observations, num_steps denote the length of the bandwidth sequence and num_alpha denotes the number of significance levels. This is a binary array and the entries are set to 1 if that observation is an outlier for that particular bandwidth and significance level.}
#' \item{\code{bw}}{The set of bandwidth values.}
#' \item{\code{gpdparas}}{The GPD parameters used. }
#' \item{\code{lookoutbw}}{The bandwidth chosen by the algorithm \code{lookout} using persistent homology.}
#'
#' @examples
#' set.seed(1)
#' x1 <- cbind.data.frame(rnorm(500), rnorm(500))
#' x2 <- cbind.data.frame(rnorm(5, mean = 10, sd = 0.2), rnorm(5, mean = 10, sd = 0.2))
#' colnames(x1) <- colnames(x2) <- c("x", "y")
#' X <- rbind.data.frame(x1, x2)
#' plot(X, pch = 19)
#' outnew <- persisting_outliers_over_alpha(X, unitize = FALSE)
#' @export

persisting_outliers_over_alpha <- function(X, alpha_min = 0.01, alpha_max = 0.1, step = 0.01, st_qq = 0.9,
                                           en_sc = sqrt(5), unitize = TRUE, num_steps = 20) {
  X <- as.data.frame(X)
  num_cols <- NCOL(X)
  num_alpha <- (alpha_max - alpha_min) / step + 1
  alpha_seq <- seq(from = alpha_min, to = alpha_max, by = step)
  if (unitize) {
    X <- unitize(X)
  }
  if (num_cols == 1L) {
    phom <- TDAstats::calculate_homology(dist(X), format = "distmat")
  } else {
    phom <- TDAstats::calculate_homology(X, dim = 0)
  }

  death_radi <- phom[, 3L]
  qq_st <- quantile(death_radi, probs = st_qq)
  qq_en <- max(death_radi) * en_sc
  bw_vals <- seq(qq_st, qq_en, length.out = num_steps)
  output <- array(0, dim = c(dim(X)[1], num_steps, num_alpha))

  q_thres <- quantile(death_radi, probs = 0.5)
  dr_thres <- death_radi[death_radi >= q_thres]
  dr_thres_diff <- diff(dr_thres)
  max_persist_ind <- which.max(dr_thres_diff)
  ind1 <- min(which(death_radi >= q_thres))
  ind <- max_persist_ind + ind1 - 1L
  bw_fixed <- death_radi[ind] * sqrt(5)
  lookoutobj1 <- lookoutliers(X, alpha = 0.05, unitize = FALSE, bw = bw_fixed)
  paras <- lookoutobj1$gpd[1:2]

  for (i in seq_along(bw_vals)) {
    lookoutobj <- lookoutliers_fixed_gpd(X, alpha = 0.05, unitize = FALSE, bw = bw_vals[i], gpd = paras)
    for (j in seq(num_alpha)) {
      outinds <- which(lookoutobj$outlier_probability < alpha_seq[j])
      output[outinds, i, j] <- 1
    }
  }

  list(
    out = output,
    bw = bw_vals,
    gpdparas = paras,
    lookoutbw = bw_fixed
  )
}
