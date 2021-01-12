
#' Plots outlier persistence for a range of significance levels.
#'
#' This function plots outlier persistence for a range of significance levels
#' using the algorithm lookout, an outlier detection method that uses
#' leave-one-out kernel density estimates and generalized Pareto distributions
#' to find outliers.
#'
#' @param object The output of the function `persisting_outliers`.
#' @param alpha The significance levels to plot.
#' @param ... Other arguments currently ignored.
#'
#' @return A ggplot object.
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
#' plot(X, pch = 19)
#' outliers <- persisting_outliers(X, unitize = FALSE)
#' autoplot(outliers)
#' @export
autoplot.persistingoutliers <- function(object, alpha = object$alpha, ...) {
  which_alpha <- (round(object$alpha, 4) %in% round(alpha, 4))
  if(all(!which_alpha))
    stop("No specified alpha values available.")
  outwts <- apply(object$out[, , which_alpha, drop = FALSE], c(1, 2), sum)
  outwtsg <- cbind.data.frame(seq(NROW(outwts)), outwts)
  colnames(outwtsg)[1] <- "Observation"
  col1 <- max(which(colSums(outwtsg) != 0))
  outwtsg <- outwtsg[, seq(col1)]

  # Long form
  dfl <- tidyr::pivot_longer(outwtsg, -Observation,
    names_to = "bw",
    values_to = "Strength"
  )
  # Add bandwidths
  dfl$Bandwidth <- object$bw[as.integer(dfl$bw)]

  # Colours
  if(length(alpha) > 1) {
    col_pal1 <- c("white", "#ffffcc", "#ffeda0", "#fed976", "#feb24c",
                  "#fd8d3c", "#fc4e2a", "#e31a1c", "#bd0026", "#800026")
  } else {
    col_pal1 <- c("white", "black")
  }

  Observation <- Bandwidth <- Strength <- NULL
  p <- ggplot2::ggplot(
      dfl,
      ggplot2::aes(x = Bandwidth, y = Observation, fill = Strength)
    ) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradientn(colours = col_pal1) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black")
    )
  if (length(alpha) == 1L) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  return(p)
}
