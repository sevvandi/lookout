#' Identifies outliers using the algorithm lookout.
#'
#' This function identifies outliers using the algorithm lookout, an outlier
#' detection method that uses leave-one-out kernel density estimates and
#' generalized Pareto distributions to find outliers.
#'
#' @details
#' The algorithm has three steps.
#'
#' \strong{Scaling.} When \code{scale = TRUE}, the data are first scaled with
#' \code{\link{mvscale}()}: each column is centred at its median, and the data
#' are rotated and scaled using the Cholesky factor of the inverse of a robust
#' (MCD) covariance estimate, so that the columns are approximately
#' uncorrelated with unit scale.
#'
#' \strong{Leave-one-out density estimates.} A kernel density estimate is
#' computed at each observation using a spherically symmetric Epanechnikov
#' kernel with support radius \code{bw}, and the contribution of the
#' observation itself is then removed to give a leave-one-out estimate. When
#' \code{bw = NULL}, the bandwidth is computed from the Euclidean minimum
#' spanning tree of the (scaled) data by \code{\link{find_tda_bw}()}: the
#' \code{gamma} quantile (type 8) of the tree's edge lengths is multiplied by
#' \code{sqrt(m + 4)}, where \code{m = NCOL(X)}. A spherically symmetric
#' Epanechnikov kernel with support radius \code{h} in \code{m} dimensions has
#' standard deviation \code{h / sqrt(m + 4)} in each coordinate, so this
#' scaling makes the kernel's per-coordinate standard deviation equal to the
#' quantile. When \code{fast = TRUE}, each kernel is summed over the \code{k}
#' nearest neighbours of the point only; see the \code{fast} argument.
#'
#' \strong{Extreme value model.} The negative logarithms of the density
#' estimates (the surprisals) above their \code{beta} quantile are modelled
#' with a generalized Pareto distribution (GPD), fitted by maximum likelihood
#' using \code{\link[evd]{fpot}()}. Because the surprisals are bounded, the
#' shape parameter is constrained to be at most zero: if the unconstrained
#' estimate is positive, the GPD is refitted with the shape fixed at zero. The
#' fitted GPD gives, for each observation, the probability of a leave-one-out
#' surprisal at least as large as the one observed, multiplied by
#' \code{1 - beta}. Observations whose probability is below \code{alpha} are
#' declared outliers.
#'
#' Setting \code{old_version = TRUE} gives the algorithm of Kandanaarachchi and
#' Hyndman (2022) instead. It differs in three places: the data are scaled so
#' that each column lies in \code{[0, 1]}, rather than with
#' \code{\link{mvscale}()}; the bandwidth is the lower end of the largest gap
#' between consecutive minimum spanning tree edge lengths among those at or
#' above their median, rather than the \code{gamma} quantile, again multiplied
#' by \code{sqrt(m + 4)}; and the GPD shape parameter is not constrained.
#'
#' @param X The numerical input data in a data.frame, matrix or tibble format.
#' @param alpha The level of significance. Default is \code{0.01}. So there is
#' a 1/100 chance of any point being falsely classified as an outlier.
#' @param beta The quantile threshold used in the GPD estimation. Default is \code{0.90}.
#' To ensure there is enough data available, values greater than 0.90 are set to 0.90.
#' @param gamma The quantile of the minimum spanning tree edge lengths used to
#' compute the bandwidth. Default is \code{0.98}. Ignored if \code{bw} is
#' provided, and ignored when \code{old_version = TRUE}, where the largest gap
#' between consecutive edge lengths is used instead. See Details.
#' @param bw The support radius of the Epanechnikov kernel, on the scale of the
#' data after any scaling. If \code{NULL} (default), it is computed from the
#' minimum spanning tree of the data as described in Details.
#' @param gpd Generalized Pareto distribution parameters. If `NULL` (the
#' default), these are estimated from the data.
#' @param scale If \code{TRUE} (the default), the data are scaled before the
#' bandwidth and density estimates are computed: with \code{\link{mvscale}()}
#' when \code{old_version = FALSE}, so that the columns are approximately
#' uncorrelated with unit scale, or by scaling each column to the range
#' \code{[0, 1]} when \code{old_version = TRUE}.
#' @param fast If \code{TRUE}, each kernel density estimate is a sum over the
#' \code{k} nearest neighbours of the point only, including the point itself,
#' where \code{k = min(max(ceiling(n / 200), 100), n, 500)} and
#' \code{n = NROW(X)}; so \code{k} is between 100 and 500, and equals 100
#' whenever \code{100 <= n <= 20000}. Wherever more than \code{k - 1} other
#' observations lie inside the kernel support, which is typical in the bulk
#' of the data in three or more dimensions, the sum is truncated and the
#' estimate is lower than the exact one. Sparse observations, which are the
#' candidates for outliers, have fewer than \code{k} neighbours inside the
#' support and their estimates are unchanged, although the GPD threshold and
#' fit can still differ. If \code{FALSE}, each kernel is summed over exactly
#' the observations inside its support, found with a fixed-radius search
#' (\code{\link[dbscan]{frNN}()}). The time and memory of this exact
#' computation are proportional to the number of pairs of observations within
#' \code{bw} of each other. In two dimensions this is usually modest, but in
#' three or more dimensions the kernel support in the bulk of the data can
#' contain thousands of observations, so the exact computation is still of
#' order \code{n^2} in the worst case and \code{fast = TRUE} remains the
#' practical choice for large \code{n}. Default is \code{TRUE} when
#' \code{NROW(X) > 10000}. The bandwidth calculation always uses all of the
#' data.
#' @param old_version If \code{TRUE}, the algorithm of Kandanaarachchi and
#' Hyndman (2022) is used. Default is \code{FALSE}, giving the algorithm of
#' Hyndman, Kandanaarachchi and Turner (2026). See Details.
#' @return A list with the following components:
#' \item{\code{data}}{The input data \code{X}, before any scaling.}
#' \item{\code{outliers}}{The set of outliers.}
#' \item{\code{outlier_probability}}{The GPD probability of the data.}
#' \item{\code{outlier_scores}}{The outlier scores of the data.}
#' \item{\code{bandwidth}}{The support radius of the Epanechnikov kernel:
#' either \code{bw}, or the value computed from the minimum spanning tree
#' multiplied by \code{sqrt(NCOL(X) + 4)}.}
#' \item{\code{kde}}{The kernel density estimate values.}
#' \item{\code{lookde}}{The leave-one-out kde values.}
#' \item{\code{gpd}}{The fitted GPD parameters.}
#' \item{\code{call}}{The matched call.}
#'@references Kandanaarachchi, S, and Hyndman, RJ (2022) Leave-one-out kernel
#' density estimates for outlier detection,
#' *J Computational & Graphical Statistics*, **31**(2), 586-599.
#' <https://robjhyndman.com/publications/lookout/>.
#'
#' Hyndman, RJ, Kandanaarachchi, S, and Turner, K (2026) Lookout 2: Anomaly
#' detection via leave-one-out kernel density estimation, arXiv:2603.22636.
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
#' @importFrom stats dist quantile median sd optimize var
#' @importFrom utils packageVersion
lookout <- function(
  X,
  alpha = 0.01,
  beta = 0.90,
  gamma = 0.98,
  bw = NULL,
  gpd = NULL,
  scale = TRUE,
  fast = NROW(X) > 10000,
  old_version = FALSE
) {
  # alpha, beta and gamma need to be between 0 and 1
  if (alpha < 0 || alpha > 1) {
    stop("alpha should be between 0 and 1.")
  }
  if (beta < 0 || beta > 1) {
    stop("beta should be between 0 and 1.")
  }
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
      gamma = gamma,
      use_differences = old_version
    ) *
      sqrt(NCOL(X) + 4)
  } else {
    bandwidth <- bw
  }

  # find kde and lookde estimates
  kdeobj <- lookde(X, bandwidth = bandwidth, fast = fast)
  log_dens <- -log(kdeobj$kde)

  # find POT GPD parameters, threshold 0.90
  beta <- min(0.9, beta)
  qq <- quantile(log_dens, probs = beta)

  # The GPD is fitted to the surprisals above their beta quantile. The largest
  # possible surprisal belongs to an observation with no other observation
  # inside its kernel support, so if more than 100(1 - beta)% of observations
  # are isolated in this way, the quantile equals the maximum and there are no
  # exceedances to fit.
  if (!any(log_dens > qq)) {
    stop(
      "Unable to fit the generalized Pareto distribution: more than ",
      sprintf("%g%%", 100 * (1 - beta)),
      " of observations have no other observation inside the kernel support ",
      "(within `bw` = ", format(bandwidth, digits = 4), " of them), so their ",
      "surprisals tie at the maximum and there are no exceedances above the ",
      "`beta` quantile. Use a larger `gamma`, or supply a larger `bw`, to ",
      "widen the kernel support."
    )
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
  m <- NCOL(x)

  # Spherically symmetric Epanechnikov kernel on the ball of radius
  # `bandwidth` in m dimensions. Its value at the origin is k0, which reduces
  # to 0.75 / bandwidth when m = 1.
  vol_unit_ball <- pi^(m / 2) / gamma(m / 2 + 1)
  k0 <- (m + 2) / (2 * vol_unit_ball * bandwidth^m)

  if (fast) {
    # Sum each kernel over the kk nearest neighbours of the point only
    # (including the point itself), where kk is between 100 and 500.
    # This truncates the sum wherever more than kk - 1 other observations lie
    # inside the kernel support, lowering the estimate in the bulk of the data.
    # Sparse observations, which have fewer than kk neighbours within the
    # support, are unaffected.
    # kNN() excludes the point itself, so ask for kk - 1 neighbours and put
    # the point back in the first column at distance 0.
    kk <- min(max(ceiling(nn / 200), 100), nn, 500)
    dist <- cbind(0, dbscan::kNN(x, k = kk - 1)$dist)
    dist[dist > bandwidth] <- NA_real_
    phat <- k0 / nn * rowSums(1 - (dist / bandwidth)^2, na.rm = TRUE)
  } else {
    # Sum each kernel over exactly the observations inside its support, found
    # with a fixed-radius search. The point itself is excluded by frNN() and
    # contributes 1 to the sum. Time and memory are proportional to the number
    # of pairs within `bandwidth` of each other rather than to nn^2.
    nbrs <- dbscan::frNN(x, eps = bandwidth, sort = FALSE)$dist
    kernel_sum <- vapply(
      nbrs,
      function(d) 1 + sum(1 - (d / bandwidth)^2),
      numeric(1L)
    )
    phat <- k0 / nn * kernel_sum
  }

  # leave one out
  kdevalsloo <- k0 / (nn - 1)
  lookde <- nn * phat / (nn - 1) - kdevalsloo

  list(x = x, kde = phat, lookde = pmax(lookde, 0))
}
