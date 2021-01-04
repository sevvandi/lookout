#' Plots outlier persistence for a given significance level.
#'
#' This function plots outlier persistence for a given significance level. using the algorithm lookout, an outlier detection method that uses leave-one-out kernel density estimates and generalized Pareto distributions to find outliers.
#'
#' @param outnew The output of the function \code{persisting_outliers_over_alpha}.
#' @param alpha The level of significance, set to \code{0.05}.
#'
#' @return A ggplot object.
#'
#' @examples
#' set.seed(1)
#' x1 <- cbind.data.frame(rnorm(500), rnorm(500))
#' x2 <- cbind.data.frame(rnorm(5, mean=10, sd=0.2), rnorm(5, mean=10, sd=0.2))
#' colnames(x1) <- colnames(x2) <- c("x", "y")
#' X <- rbind.data.frame(x1, x2)
#' plot(X, pch=19)
#' outnew <- persisting_outliers_over_alpha(X, unitize=FALSE)
#' plot_outlier_barcode(outnew, alpha=0.05)
#'
#' @export plot_outlier_barcode
plot_outlier_barcode <- function(outnew, alpha=0.05){
  ind <- alpha*100
  mat <-  outnew$out[ , ,ind]
  mat2 <- cbind(1:dim(mat)[1], mat)
  mat2 <- as.data.frame(mat2)
  colnames(mat2)[1] <- "id"
  colnames(mat2)[2:(length(outnew$bw)+1)] <- outnew$bw

  colsums <- apply(mat2, 2, sum)
  col1 <- max(which(colsums!=0)) + 1
  col2 <- dim(mat2)[2]
  mat2 <- mat2[ ,-c(col1:col2)]

  dfl <- tidyr::pivot_longer(mat2, 2:dim(mat2)[2])
  colnames(dfl) <- c("x", "y", "outlier")
  dfl$y <- as.numeric(paste(dfl$y))
  ggplot2::ggplot(dfl,ggplot2::aes(x=dfl$y,y=dfl$x, fill=dfl$outlier)) + ggplot2::geom_raster() + ggplot2::scale_fill_gradientn(colours=c("white", "black")) + ggplot2::ylab("Observations") + ggplot2::xlab("Bandwidth")+  ggplot2::theme(legend.position = "none") + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))
}

#' Plots outlier persistence for a range of significance levels.
#'
#' This function plots outlier persistence a range of significance levels.

#' @inheritParams plot_outlier_barcode
#'
#' @return A ggplot object.
#'
#' @examples
#' set.seed(1)
#' x1 <- cbind.data.frame(rnorm(500), rnorm(500))
#' x2 <- cbind.data.frame(rnorm(5, mean=10, sd=0.2), rnorm(5, mean=10, sd=0.2))
#' colnames(x1) <- colnames(x2) <- c("x", "y")
#' X <- rbind.data.frame(x1, x2)
#' plot(X, pch=19)
#' outnew <- persisting_outliers_over_alpha(X, unitize=FALSE)
#' plot_outlier_barcode_over_alpha(outnew)
#'
#' @export
plot_outlier_barcode_over_alpha <- function(outnew){
  col_pal1 <- c("white", "#ffffcc", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#bd0026", "#800026")
  outwts <- apply(outnew$out, c(1,2), sum)
  outwtsg <- cbind.data.frame(1:dim(outwts)[1], outwts)
  colnames(outwtsg)[1] <- "id"
  colnames(outwtsg)[2:21] <- outnew$bw
  colsums <- apply(outwtsg, 2, sum)
  col1 <- max(which(colsums!=0)) + 1
  col2 <- dim(outwtsg)[2]
  outwtsg <- outwtsg[ ,-c(col1:col2)]


  dfl <- tidyr::pivot_longer(outwtsg, 2:dim(outwtsg)[2])
  colnames(dfl) <- c("x", "y", "strength")
  dfl$y <- as.numeric(paste(dfl$y))
  ggplot2::ggplot(dfl,ggplot2::aes(x=dfl$y,y=dfl$x, fill=dfl$strength)) + ggplot2::geom_raster() + ggplot2::scale_fill_gradientn(colours=col_pal1)+ ggplot2::ylab("Observations") + ggplot2::xlab("Bandwidth") + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))
}
