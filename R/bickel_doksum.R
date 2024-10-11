## Bickel-Doksum transformations

bickeldoksum <- function(x, ...) {
  stopifnot(is.numeric(x))
  l <- estimate_bickeldoksum_lambda(x, ...)
  bickeldoksum_trans(x,l)
}

estimate_bickeldoksum_lambda <- function (x, lower = 0, upper = 1, eps = 0.00001)
{
  n <- length(x)
  x <- x[!is.na(x)]
  bd_loglik <- function(lambda) {
    z <- bickeldoksum_trans(x, lambda, eps)
    var_z <- var(z) * (n - 1)/n
    -0.5 * n * log(var_z)
  }
  results <- optimize(bd_loglik, lower = lower, upper = upper,
                      maximum = TRUE, tol = 1e-04)
  results$maximum
}

bickeldoksum_trans <- function (x, lambda, eps = 0.001)
{
  if (lambda < 0)
    x[x < 0] <- NA
  if (abs(lambda) < eps)
    val <- log(x)
  else val <- (sign(x) * abs(x)^lambda - 1)/lambda
  val
}


