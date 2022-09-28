#' Identifies outliers using the algorithm lookout.
#'
#' This function identifies outliers using the algorithm lookout, an outlier
#' detection method that uses leave-one-out kernel density estimates and
#' generalized Pareto distributions to find outliers.
#'
#' @param X The input data in a dataframe, matrix or tibble format.
#' @param alpha The level of significance. Default is \code{0.05}.
#' @param unitize An option to normalize the data. Default is \code{TRUE},
#'   which normalizes each column to \code{[0,1]}.
#' @param bw Bandwidth parameter. Default is \code{NULL} as the bandwidth is
#'   found using Persistent Homology.
#' @param gpd Generalized Pareto distribution parameters. If `NULL` (the
#' default), these are estimated from the data.
#' @param fast If set to \code{TRUE}, makes the computation faster by sub-setting
#'   the data for the bandwidth calculation.
#'
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
#'   data.frame(x = rnorm(500),
#'              y = rnorm(500)),
#'   data.frame(x = rnorm(5, mean = 10, sd = 0.2),
#'              y = rnorm(5, mean = 10, sd = 0.2))
#' )
#' lo <- lookout(X)
#' lo
#' autoplot(lo)
#' @export lookout
#' @importFrom stats dist quantile median
lookout <- function(X,
                    alpha = 0.05,
                    unitize = TRUE,
                    bw = NULL,
                    gpd = NULL,
                    fast = TRUE) {
  # Prepare X matrix
  origX <- X
  X <- as.matrix(X)
  if (unitize) {
    X <- unitize(X)
  }

  # Find bandwidth and scale for Epanechnikov kernel
  if (is.null(bw)) {
    bandwidth <- find_tda_bw(X, fast = fast) * sqrt(5)
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

  if(is.null(gpd)) {
    M1 <- evd::fpot(log_dens, qq, std.err = FALSE)
    gpd <- M1$estimate[1L:2L]
  }
  # for these Generalized Pareto distribution parameters, compute the
  # probabilities of leave-one-out kernel density estimates
  potlookde <- evd::pgpd(-log(kdeobj$lookde), loc = qq,
    scale = gpd[1], shape = gpd[2], lower.tail = FALSE)
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
  ), class='lookoutliers')
}

find_tda_bw <- function(X, fast) {
  X <- as.matrix(X)

  # select a subset of X for tda computation
  if(fast){
    Xsub <- subset_for_tda(X)
  }else{
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
  dr_thres_diff <- diff(death_radi_upper)
  return(death_radi_upper[which.max(dr_thres_diff)])
}

lookde <- function(x, bandwidth, fast) {
  x <- as.matrix(x)
  nn <- NROW(x)

  if(fast){
    # To make the nearest neighbour distance computation faster
    # select a kk different to nn as follows
    kk <- min(max(ceiling(nn/200), 100), nn, 500)
  }else{
    kk <- nn
  }

  # Epanechnikov kernel density estimate
  dist <- RANN::nn2(x, k = kk)$nn.dists
  dist[dist > bandwidth] <- NA_real_
  phat <- 0.75 / (nn*bandwidth) * rowSums(1-(dist/bandwidth)^2, na.rm=TRUE)

  # leave one out
  kdevalsloo <- 0.75 / ((nn - 1) * (bandwidth))
  lookde <- nn * phat / (nn - 1) - kdevalsloo

  list(x = x, kde = phat, lookde = pmax(lookde, 0))
}


subset_for_tda <- function(X){
  # Leader algorithm in HDoutliers
  # Inserted from HDoutliers function getHDmembers
  # We cannot call that function because the algorithm only comes to
  # effect if the number of rows are greater than 10000
  # And we have used RANN::nn2, which is a faster algorithm.

  X <- as.matrix(X)

  n <- nrow(X)
  p <- ncol(X)
  radius <- 0.1/(log(n)^(1/p))
  members <- rep(list(NULL), n)
  exemplars <- 1
  members[[1]] <- 1

  for (i in 2:n) {
    KNN <- RANN::nn2(data = X[c(exemplars,i), , drop = F],
                     query = X[i, , drop = F], k = 2)
    m <- KNN$nn.idx[1, 2]
    d <- KNN$nn.dists[1, 2]
    if (d < radius) {
      l <- exemplars[m]
      members[[l]] <- c(members[[l]], i)
      next
    }
    exemplars <- c(exemplars, i)
    members[[i]] <- i
  }
  X[exemplars, ]
}
