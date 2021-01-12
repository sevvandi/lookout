#' @method print persistingoutliers
#' @export

print.persistingoutliers <- function(x, ...) {

  cat("Persistent outliers using lookout algorithm")
  cat("\n\nCall: ")
  print(x$call)
  cat("\nLookout bandwidth: ", x$lookoutbw,"\n")
}

#' @method print lookoutliers
#' @export
print.lookoutliers <- function(x, ...) {
  cat("Leave-out-out KDE outliers using lookout algorithm")
  cat("\n\nCall: ")
  print(x$call)
  cat("\n")
  print(x$outliers)
  cat("\n")
}
