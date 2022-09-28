#' Plots outliers identified by lookout algorithm.
#'
#' Scatterplot of two columns from the data set with outliers highlighted.
#'
#' @param object The output of the function `lookout`.
#' @param columns Which columns of the original data to plot
#'  (specified as either numbers or strings)
#' @param ... Other arguments currently ignored.
#'
#' @return A ggplot object.
#'
#' @examples
#' X <- rbind(
#'   data.frame(x = rnorm(500),
#'              y = rnorm(500)),
#'   data.frame(x = rnorm(5, mean = 10, sd = 0.2),
#'              y = rnorm(5, mean = 10, sd = 0.2))
#' )
#' lo <- lookout(X)
#' autoplot(lo)
#' @export
autoplot.lookoutliers <- function(object, columns=1:2, ...) {
  # Column names
  varnames <- colnames(object$data)
  if(is.null(varnames))
    varnames <- paste0("V",seq(NCOL(object$data)))
  X <- as.data.frame(object$data)
  colnames(X) <- varnames
  if(is.character(columns)) {
    columns <- match(columns, varnames)
  } else {
    columns <- columns[columns <= NCOL(X)]
  }

  # Outliers
  outliers <- NULL
  X$outliers <- rep(FALSE, NROW(X))
  X$outliers[object$outliers[,"Outliers"]] <- TRUE

  # y axis
  if(length(columns) > 1) {
    ..y <- X[,columns[2L]]
    ..yvar <- varnames[columns[2L]]
  } else {
    ..y <- 0
    ..yvar <- ""
  }

  # Produce plot
  p <- ggplot2::ggplot(X, ggplot2::aes(x=X[,columns[1L]], y=..y)) +
      ggplot2::geom_point(ggplot2::aes(col=outliers)) +
      ggplot2::labs(x=varnames[columns[1L]], y=..yvar) +
      ggplot2::scale_color_manual(values=c(`TRUE` = "red", `FALSE`="black")) +
      ggplot2::guides(color = "none")
  if(NCOL(object$data) == 1L) {
    p <- p + ggplot2::scale_y_continuous(breaks=NULL)
  }
  p
}
