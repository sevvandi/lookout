#' Identifies outliers using the algorithm lookout.
#'
#' This function identifies outliers using the algorithm lookout, an outlier detection method that uses leave-one-out kernel density estimates and generalized Pareto distributions to find outliers.
#'
#' @param X The input data in a dataframe, matrix or tibble format.
#' @param alpha The level of significance. Default is \code{0.0.05}.
#' @param unitize An option to normalize the data. Default is \code{TRUE}, which normalizes each column to \code{[0,1]}.
#' @param bw Bandwidth parameter. Default is \code{NULL} as the bandwidth is found using Persistent Homology.
#'
#' @return A list with the following components:
#' \item{\code{outliers}}{The set of outliers.}
#' \item{\code{outlier_probability}}{The GPD probability of the data.}
#' \item{\code{bandwidth}}{The bandwdith selected using persistent homology. }
#' \item{\code{kde}}{The kernel density estimate values.}
#' \item{\code{lookde}}{The leave-one-out kde values.}
#' \item{\code{gpd}}{The fitted GPD parameters.}
#'
#' @examples
#' set.seed(1)
#' x1 <- cbind.data.frame(rnorm(500), rnorm(500))
#' x2 <- cbind.data.frame(rnorm(5, mean=10, sd=0.2), rnorm(5, mean=10, sd=0.2))
#' colnames(x1) <- colnames(x2) <- c("x", "y")
#' X <- rbind.data.frame(x1, x2)
#' plot(X, pch=19)
#' out <- lookoutliers(X)
#'
#' @export lookoutliers
#' @importFrom stats dist quantile
lookoutliers <- function(X, alpha=0.05,  unitize=TRUE, bw=NULL){
  X <- as.data.frame(X)
  dd <- dim(X)[2]
  if(unitize){
    for(col in 1:dd){
      X[ ,col] <- (X[ ,col] - min(X[ ,col]) )/(max(X[ ,col]) - min(X[ ,col]))
    }
  }
  if(is.null(bw)){
    # find bandwidth
    bandwidth <- find_tda_bw(X)
    # scale it for Epanechnikov kernel
    bandwidth <- sqrt(5)*bandwidth
  }else{
    bandwidth <- bw
  }

  # find kde and lookde estimates
  kdeobj <- lookde(X, bandwidth = bandwidth)
  kdevals <- kdeobj$kde
  lookde <- kdeobj$lookde
  # take logs
  log_dens <- -1*log(kdevals)
  # find POT GPD parameters, threshold 0.90
  qq <- quantile(log_dens, probs=min(0.9, 1-alpha))
  # check if there are points above the quantile
  len_over <- length(which(log_dens > qq))
  if(len_over==0){
    rpts <- 1
    while(len_over==0){
      qq <- quantile(log_dens, probs=(0.9 - rpts*0.05))
      len_over <- length(which(log_dens > qq))
      rpts <- rpts + 1
    }
  }
  M1 <- evd::fpot(log_dens, qq, std.err = FALSE)
  # for these Generalized Pareto distribution paramaters, compute the probabilities of leave-one-out kernel density estimates
  loglookde <- -1*log(lookde)
  potlookde <- evd::pgpd(loglookde, loc=qq,  scale= M1$estimate[1], shape=M1$estimate[2], lower.tail = FALSE )
  # select outliers according to threshold
  outliers <- which(potlookde < alpha)
  dfout <- cbind.data.frame(outliers, potlookde[outliers])
  colnames(dfout) <- c("Outliers", "Probability")

  out <- list()


  out$outliers <- dfout
  out$outlier_probability <- potlookde
  out$bandwidth <- bandwidth
  out$kde <- kdevals
  out$lookde <- lookde
  out$gpd <- c(M1$estimate[1], M1$estimate[2], M1$estimate[3])

  return(out)
}

find_tda_bw <- function(X){
  X <- as.data.frame(X)
  num_cols <- dim(X)[2]
  if(num_cols==1){
    distx <- dist(X)
    phom <- TDAstats::calculate_homology(distx, format="distmat")
  }else{
    phom <- TDAstats::calculate_homology(X, dim=0)
  }

  death_radi <- phom[ ,3]
  # q_thres <- quantile(death_radi, probs=0)
  # dr_thres <- death_radi[death_radi >= q_thres]
  # dr_thres_diff <- diff(dr_thres)
  dr_thres_diff <- diff(death_radi)
  max_persist_ind <- which.max(dr_thres_diff)
  #ind1 <- min(which(death_radi >= q_thres))
  # ind <- max_persist_ind + ind1 -1
  ind <- max_persist_ind #-1
  bw <- death_radi[ind]
  return(bw)
}


lookde <- function(x, bandwidth){
  x <- as.data.frame(x)
  dd <- dim(x)[2]
  nn <- dim(x)[1]
  k_val <- max(100, ceiling(nn/3))
  nnobj <- RANN::nn2(x, k=nn)  # k_val
  phat <- rep(0, nn)
  lookde <- rep(0, nn)
  for(kk in 1:nn){
    inds <- which(nnobj$nn.dists[kk, ] < bandwidth)
    kdevals <- 0.75/(nn*(bandwidth))*(1-nnobj$nn.dists[kk, inds]^2/bandwidth^2)
    phat[kk] <- sum(kdevals)
  }

  # leave one out
  kdevalsloo <- 0.75/((nn-1)*(bandwidth))
  lookde <- nn*phat/(nn-1) - kdevalsloo

  neginds <- which(lookde < 0)
  lookde[neginds] <- 0

  out <- list()
  out$x <- x
  out$kde <- phat
  out$lookde <- lookde
  return(out)
}


lookoutliers_fixed_gpd <- function(X, alpha=0.05, unitize=TRUE, bw=NULL, gpd){
  X <- as.data.frame(X)
  dd <- dim(X)[2]
  if(unitize){
    for(col in 1:dd){
      X[ ,col] <- (X[ ,col] - min(X[ ,col]) )/(max(X[ ,col]) - min(X[ ,col]))
    }
  }
  if(is.null(bw)){
    # find bandwidth
    bandwidth <- find_tda_bw(X)
    # scale it for Epanechnikov kernel
    bandwidth <- sqrt(5)*bandwidth
  }else{
    bandwidth <- bw
  }

  # find kde and lookde estimates
  kdeobj <- lookde(X, bandwidth = bandwidth)
  kdevals <- kdeobj$kde
  lookde <- kdeobj$lookde
  # take logs
  log_dens <- -1*log(kdevals)
  # find POT GPD parameters, threshold 0.90
  qq <- quantile(log_dens, probs=min(0.9, 1-alpha))
  # check if there are points above the quantile
  len_over <- length(which(log_dens > qq))
  if(len_over==0){
    rpts <- 1
    while(len_over==0){
      qq <- quantile(log_dens, probs=(0.9 - rpts*0.05))
      len_over <- length(which(log_dens > qq))
      rpts <- rpts + 1
    }
  }
  # M1 <- evd::fpot(log_dens, qq, std.err = FALSE)
  # M1 <- evd::fpot(log_dens, qq, scale=gpd[1], shape= gpd[2], std.err = FALSE)
  # for these Generalized Pareto distribution paramaters, compute the probabilities of leave-one-out kernel density estimates
  loglookde <- -1*log(lookde)
  # potlookde <- evd::pgpd(loglookde, loc=qq,  scale= M1$estimate[1], shape=gpd, lower.tail = FALSE )
  potlookde <- evd::pgpd(loglookde, loc=qq,  scale=gpd[1], shape=gpd[2], lower.tail = FALSE )
  # select outliers according to threshold
  outliers <- which(potlookde < alpha)
  dfout <- cbind.data.frame(outliers, potlookde[outliers])
  colnames(dfout) <- c("Outliers", "Probability")

  out <- list()


  out$outliers <- dfout
  out$outlier_probability <- potlookde
  out$bandwidth <- bandwidth
  out$kde <- kdevals
  out$lookde <- lookde
  out$gpd <- c(gpd[1], gpd[2])

  return(out)
}
