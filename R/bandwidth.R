#' Computes the bandwidth for lookout from the minimum spanning tree.
#'
#' This function computes the bandwidth used by \code{\link{lookout}()} for
#' its kernel density estimates, from the edge lengths of the Euclidean
#' minimum spanning tree of the data.
#'
#' @details The Euclidean minimum spanning tree of the rows of \code{X} is
#' computed with \code{\link[mlpack]{emst}()}, giving \code{NROW(X) - 1} edge
#' lengths in increasing order. Each edge length is the distance at which two
#' clusters merge under single-linkage clustering, or equivalently the death
#' radius of a connected component in the Rips filtration of the data, which
#' is why the function is named after topological data analysis.
#'
#' By default (\code{use_differences = FALSE}), the value returned is the
#' \code{gamma} quantile of the edge lengths, computed with
#' \code{\link[stats]{quantile}()} using \code{type = 8}. With
#' \code{use_differences = TRUE}, only the edge lengths at or above their
#' median are kept, and the value returned is the lower of the two consecutive
#' edge lengths with the largest gap between them. This is the rule of
#' Kandanaarachchi and Hyndman (2022), used by \code{lookout()} when
#' \code{old_version = TRUE}.
#'
#' The value returned is on the scale of \code{X} and is not itself the
#' support radius of the kernel. \code{\link{lookout}()} multiplies it by
#' \code{sqrt(m + 4)}, where \code{m = NCOL(X)}, to obtain the support radius
#' of the Epanechnikov kernel. A spherically symmetric Epanechnikov kernel with
#' support radius \code{h} in \code{m} dimensions has standard deviation
#' \code{h / sqrt(m + 4)} in each coordinate, so the kernel used by
#' \code{lookout()} has per-coordinate standard deviation equal to the value
#' returned here. Note that \code{lookout()} applies this function to the
#' scaled data when \code{scale = TRUE}.
#'
#' The minimum spanning tree is always computed on all of \code{X}.
#'
#' @inheritParams lookout
#' @param fast Deprecated and ignored. The bandwidth is always computed
#' using all of \code{X}.
#' @param gamma The quantile of the minimum spanning tree edge lengths to
#' return. Default is \code{0.98}. Ignored when
#' \code{use_differences = TRUE}.
#' @param use_differences If \code{TRUE}, the value returned is the lower end
#' of the largest gap between consecutive minimum spanning tree edge lengths
#' among those at or above their median. If \code{FALSE} (the default), the
#' \code{gamma} quantile of the edge lengths is returned.
#'
#' @return A single number: the bandwidth on the scale of \code{X}, before
#' multiplication by \code{sqrt(NCOL(X) + 4)}.
#'
#' @references Kandanaarachchi, S, and Hyndman, RJ (2022) Leave-one-out kernel
#' density estimates for outlier detection,
#' *J Computational & Graphical Statistics*, **31**(2), 586-599.
#' <https://robjhyndman.com/publications/lookout/>.
#'
#' Hyndman, RJ, Kandanaarachchi, S, and Turner, K (2026) Lookout 2: Anomaly
#' detection via leave-one-out kernel density estimation, arXiv:2603.22636.
#' <https://robjhyndman.com/publications/lookout2.html>
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
