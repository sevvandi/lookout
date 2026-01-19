#' Identifies outliers using the algorithm lookout.
#'
#' This function identifies outliers using the algorithm lookout, an outlier
#' detection method that uses leave-one-out kernel density estimates and
#' generalized Pareto distributions to find outliers.
#'
#' @param X The numerical input data in a data.frame, matrix or tibble format.
#' @param alpha The level of significance. Default is \code{0.01}. So there is
#' a 1/100 chance of any point being falsely classified as an outlier.
#' @param beta The quantile threshold used in the GPD estimation. Default is \code{0.90}.
#' To ensure there is enough data available, values greater than 0.90 are set to 0.90.
#' @param gamma Parameter for bandwidth calculation giving the quantile of the
#' Rips death radii to use for the bandwidth. Default is \code{0.97}. Ignored
#' under the old version; where the lower limit of the maximum Rips death radii
#' difference is used. Also ignored if \code{bw} is provided.
#' @param bw Bandwidth parameter. If \code{NULL} (default), the bandwidth is
#'   found using Persistent Homology.
#' @param gpd Generalized Pareto distribution parameters. If `NULL` (the
#' default), these are estimated from the data.
#' @param scale If \code{TRUE}, the data is standardized. Using the old version,
#' unit scaling is applied so that each column is in the range \code{[0,1]}.
#' Under the new version, robust rotation and scaling is used so that the columns
#' are approximately uncorrelated with unit variance. Default is \code{TRUE}.
#' @param fast If \code{TRUE} (default), makes the computation faster by
#' sub-setting the data for the bandwidth calculation.
#' @param old_version Logical indicator of which version of the algorithm to use.
#' Default is FALSE, meaning the newer version is used.
#' @return A list with the following components:
#' \item{\code{outliers}}{The set of outliers.}
#' \item{\code{outlier_probability}}{The GPD probability of the data.}
#' \item{\code{outlier_scores}}{The outlier scores of the data.}
#' \item{\code{bandwidth}}{The bandwdith selected using persistent homology. }
#' \item{\code{kde}}{The kernel density estimate values.}
#' \item{\code{lookde}}{The leave-one-out kde values.}
#' \item{\code{gpd}}{The fitted GPD parameters.}
#'@references Kandanaarachchi, S, and Hyndman, RJ (2022) Leave-one-out kernel
#' density estimates for outlier detection,
#' *J Computational & Graphical Statistics*, **31**(2), 586-599.
#' <https://robjhyndman.com/publications/lookout/>.
#'
#' Hyndman, RJ, Kandanaarachchi, S, and Turner, K (2026) When lookout meets
#' crackle: Anomaly detection using kernel density estimation, in preparation.
#' <https://robjhyndman.com/publications/lookout2.html>
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
#' lo <- lookout(X)
#' lo
#' autoplot(lo)
#' @export lookout
#' @importFrom stats dist quantile median sd
lookout <- function(
  X,
  alpha = 0.01,
  beta = 0.90,
  gamma = 0.97,
  bw = NULL,
  gpd = NULL,
  scale = TRUE,
  fast = NROW(X) > 1000,
  old_version = FALSE
) {
  # alpha, beta and gamma need to be between 0 and 1
  if (alpha < 0 || alpha > 1) {
    stop("gamma should be between 0 and 1.")
  }
  if (beta < 0 || beta > 1) {
    stop("gamma should be between 0 and 1.")
  }
  # gamma needs to be between 0 and 1
  if (gamma < 0 || gamma > 1) {
    stop("gamma should be between 0 and 1.")
  }

  # Prepare X matrix
  origX <- X
  X <- as.matrix(X)
  if (scale) {
    if (old_version) {
      X <- unitize(X)
    } else {
      X <- mvscale(X)
    }
  }

  # Find bandwidth and scale for Epanechnikov kernel
  if (is.null(bw)) {
    bandwidth <- find_tda_bw(
      X,
      fast = fast,
      gamma,
      use_differences = old_version
    ) *
      sqrt(5)
  } else {
    bandwidth <- bw
  }

  # find kde and lookde estimates
  kdeobj <- lookde(X, bandwidth = bandwidth, fast = fast)
  log_dens <- -log(kdeobj$kde)

  # find POT GPD parameters, threshold 0.90
  beta <- min(0.9, beta)
  qq <- quantile(log_dens, probs = beta)

  # check if there are points above the quantile
  if (!any(log_dens > qq)) {
    stop("No points above the quantile for GPD estimation")
  }

  if (is.null(gpd)) {
    M1 <- evd::fpot(log_dens, qq, std.err = FALSE)
    gpd <- M1$estimate[1L:2L]
    if (gpd[2] > 0 & !old_version) {
      # This should only be done in the new lookout
      # This shows that shape is estimated to be positive.
      # This should not be the case because log densities are bounded
      M1 <- evd::fpot(log_dens, qq, shape = 0, std.err = FALSE)
      gpd <- c(M1$estimate, 0)
    }
  }
  # for these Generalized Pareto distribution parameters, compute the
  # probabilities of leave-one-out kernel density estimates
  potlookde <- evd::pgpd(
    -log(kdeobj$lookde),
    loc = qq,
    scale = gpd[1],
    shape = gpd[2],
    lower.tail = FALSE
  ) *
    (1 - beta)

  outscores <- 1 - potlookde
  # select outliers according to threshold
  outliers <- which(potlookde < alpha)
  dfout <- cbind.data.frame(outliers, potlookde[outliers])
  colnames(dfout) <- c("Outliers", "Probability")

  structure(
    list(
      data = origX,
      outliers = dfout,
      outlier_probability = potlookde,
      outlier_scores = outscores,
      bandwidth = bandwidth,
      kde = kdeobj$kde,
      lookde = kdeobj$lookde,
      gpd = gpd,
      call = match.call()
    ),
    class = "lookoutliers"
  )
}


lookde <- function(x, bandwidth, fast) {
  x <- as.matrix(x)
  nn <- NROW(x)

  if (fast) {
    # To make the nearest neighbour distance computation faster
    # select a kk different to nn as follows
    kk <- min(max(ceiling(nn / 200), 100), nn, 500)
  } else {
    kk <- nn
  }

  # Epanechnikov kernel density estimate
  dist <- RANN::nn2(x, k = kk)$nn.dists
  dist[dist > bandwidth] <- NA_real_
  phat <- 0.75 /
    (nn * bandwidth) *
    rowSums(1 - (dist / bandwidth)^2, na.rm = TRUE)

  # leave one out
  kdevalsloo <- 0.75 / ((nn - 1) * (bandwidth))
  lookde <- nn * phat / (nn - 1) - kdevalsloo

  list(x = x, kde = phat, lookde = pmax(lookde, 0))
}


subset_for_tda <- function(X) {
  # Leader algorithm in HDoutliers
  # Inserted from HDoutliers function getHDmembers
  # We cannot call that function because the algorithm only comes to
  # effect if the number of rows are greater than 10000
  # And we have used RANN::nn2, which is a faster algorithm.

  X <- as.matrix(X)

  n <- nrow(X)
  p <- ncol(X)

  Xu <- unitize(X)

  sds <- apply(Xu, 2, sd)
  sd_radius <- sqrt(sum(sds^2))
  radius <- min(0.1 / (log(n)^(1 / p)), sd_radius)
  members <- rep(list(NULL), n)
  exemplars <- 1
  members[[1]] <- 1

  for (i in 2:n) {
    KNN <- RANN::nn2(
      data = Xu[c(exemplars, i), , drop = FALSE],
      query = Xu[i, , drop = FALSE],
      k = 2
    )
    m <- KNN$nn.idx[1, 2]
    d <- KNN$nn.dists[1, 2]
    if (d < radius) {
      curr <- length(exemplars)
      l <- exemplars[curr]
      members[[l]] <- c(members[[l]], i)
      next
    }
    exemplars <- c(exemplars, i)
    members[[i]] <- i
  }
  # X[exemplars, ]
  exemplars
}
