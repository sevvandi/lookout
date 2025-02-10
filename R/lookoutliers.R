#' Identifies outliers using the algorithm lookout.
#'
#' This function identifies outliers using the algorithm lookout, an outlier
#' detection method that uses leave-one-out kernel density estimates and
#' generalized Pareto distributions to find outliers.
#'
#' @param X The numerical input data in a data.frame, matrix or tibble format.
#' @param alpha The level of significance. Default is \code{0.05}.
#' @param unitize If \code{TRUE}, the data is standardized so that
#' each column is in the range \code{[0,1]}. Default is \code{FALSE}.
#' @param normalize If \code{TRUE} (default), each column of the data is
#' transformed to be closer to a standard normal distribution, and the data
#' is rotated so that the columns are pairwise uncorrelated. This is done
#' robustly using \code{\link{weird}{mvscale}}.
#' @param bw Bandwidth parameter. If \code{NULL} (default), the bandwidth is
#'   found using Persistent Homology.
#' @param gpd Generalized Pareto distribution parameters. If `NULL` (the
#' default), these are estimated from the data.
#' @param fast If \code{TRUE} (default), makes the computation faster by
#' sub-setting the data for the bandwidth calculation.
#' @param bw_para Parameter for bandwidth calculation. Default is \code{0.95}.
#'   If set to 1, then the bandwidth corresponds to the maximum Rips death radii
#'   difference. If set to 0.95, then the bandwidth corresponds to the 95th
#'   quantile of Rips death radii. Other probabilities can be used.
#' @param transformation Ignored if \code{normalize = FALSE}. Specifies
#'  either \code{YJ} for a Yeo-Johnson transformation, or \code{BD} for a
#'  Bickel-Doksum transformation
#' @return A list with the following components:
#' \item{\code{outliers}}{The set of outliers.}
#' \item{\code{outlier_probability}}{The GPD probability of the data.}
#' \item{\code{outlier_scores}}{The outlier scores of the data.}
#' \item{\code{bandwidth}}{The bandwdith selected using persistent homology. }
#' \item{\code{kde}}{The kernel density estimate values.}
#' \item{\code{lookde}}{The leave-one-out kde values.}
#' \item{\code{gpd}}{The fitted GPD parameters.}
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
#' lo <- lookout(X)
#' lo
#' autoplot(lo)
#' @export lookout
#' @importFrom stats dist quantile median sd
lookout <- function(X,
                    alpha = 0.05,
                    unitize = TRUE,
                    normalize = FALSE,
                    bw = NULL,
                    gpd = NULL,
                    fast = NROW(X)>1000,
                    bw_para = 0.98,
                    transformation = c("YJ","BD")) {
  transformation <- match.arg(transformation)

  # bw_para needs to be between 0 and 1
  if (bw_para < 0 || bw_para > 1) {
    stop("bw_para should be between 0 and 1.")
  }

  # Prepare X matrix
  origX <- X
  X <- as.matrix(X)
  if (unitize) {
    X <- unitize(X)
  }
  if (normalize) {
    X <- transform_normal(X, transformation = transformation)
  }

  # Find bandwidth and scale for Epanechnikov kernel
  if (is.null(bw)) {
    bandwidth <- find_tda_bw(X, fast = fast, bw_para) * sqrt(5)
  } else {
    bandwidth <- bw
  }

  # find kde and lookde estimates
  kdeobj <- lookde(X, bandwidth = bandwidth, fast = fast)
  log_dens <- -log(kdeobj$kde)

  # find POT GPD parameters, threshold 0.90
  qq <- quantile(log_dens, probs = min(0.9, 1 - alpha))

  # check if there are points above the quantile
  len_over <- length(which(log_dens > qq))
  if (len_over == 0L) {
    rpts <- 1L
    while (len_over == 0L) {
      qq <- quantile(log_dens, probs = (0.9 - rpts * 0.05))
      len_over <- length(which(log_dens > qq))
      rpts <- rpts + 1L
    }
  }


  if (is.null(gpd)) {
    M1 <- evd::fpot(log_dens, qq, std.err = FALSE)
    gpd <- M1$estimate[1L:2L]
    if(gpd[2] > 0){
      # This shows that shape is estimated to be positive.
      # This should not be the case because log densities are bounded
      M1 <- evd::fpot(log_dens, qq, shape = 0, std.err = FALSE)
      gpd <- c(M1$estimate, 0)
    }
  }
  # for these Generalized Pareto distribution parameters, compute the
  # probabilities of leave-one-out kernel density estimates
  potlookde <- evd::pgpd(-log(kdeobj$lookde),
                         loc = qq,
                         scale = gpd[1], shape = gpd[2], lower.tail = FALSE
  )

  outscores <- 1 - potlookde
  # select outliers according to threshold
  outliers <- which(potlookde < alpha)
  dfout <- cbind.data.frame(outliers, potlookde[outliers])
  colnames(dfout) <- c("Outliers", "Probability")

  structure(list(
    data = origX,
    outliers = dfout,
    outlier_probability = potlookde,
    outlier_scores = outscores,
    bandwidth = bandwidth,
    kde = kdeobj$kde,
    lookde = kdeobj$lookde,
    gpd = gpd,
    call = match.call()
  ), class = "lookoutliers")
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
  phat <- 0.75 / (nn*bandwidth) * rowSums(1 - (dist/bandwidth)^2, na.rm = TRUE)

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
      query = Xu[i, , drop = FALSE], k = 2
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
  #X[exemplars, ]
  exemplars
}
