set.seed(42)
X_clean <- data.frame(x = rnorm(200), y = rnorm(200))
X_outliers <- rbind(
  X_clean,
  data.frame(x = c(10, -10), y = c(10, -10))
)

# --- lookout() ---

test_that("lookout returns correct structure", {
  lo <- lookout(X_clean)
  expect_s3_class(lo, "lookoutliers")
  expect_named(
    lo,
    c(
      "data",
      "outliers",
      "outlier_probability",
      "outlier_scores",
      "bandwidth",
      "kde",
      "lookde",
      "gpd",
      "call"
    )
  )
})

test_that("lookout detects obvious outliers", {
  lo <- lookout(X_outliers)
  n <- nrow(X_outliers)
  # The two injected outliers (rows 201, 202) should be flagged
  expect_true(201 %in% lo$outliers$Outliers)
  expect_true(202 %in% lo$outliers$Outliers)
})

test_that("lookout output lengths are consistent with input", {
  lo <- lookout(X_clean)
  n <- nrow(X_clean)
  expect_length(lo$kde, n)
  expect_length(lo$lookde, n)
  expect_length(lo$outlier_probability, n)
  expect_length(lo$outlier_scores, n)
})

test_that("lookout probabilities are in [0, 1]", {
  lo <- lookout(X_clean)
  expect_true(all(lo$outlier_probability >= 0))
  expect_true(all(lo$outlier_probability <= 1))
})

test_that("lookout scores are in [0, 1]", {
  lo <- lookout(X_clean)
  expect_true(all(lo$outlier_scores >= 0))
  expect_true(all(lo$outlier_scores <= 1))
})

test_that("lookout kde values are non-negative", {
  lo <- lookout(X_clean)
  expect_true(all(lo$kde >= 0))
  expect_true(all(lo$lookde >= 0))
})

test_that("lookout accepts matrix input", {
  lo <- lookout(as.matrix(X_clean))
  expect_s3_class(lo, "lookoutliers")
})

test_that("lookout bandwidth can be supplied manually", {
  lo <- lookout(X_clean, bw = 0.5)
  expect_equal(lo$bandwidth, 0.5)
})

test_that("lookout respects alpha threshold", {
  lo_strict <- lookout(X_outliers, alpha = 1e-6)
  lo_loose <- lookout(X_outliers, alpha = 0.5)
  expect_lte(nrow(lo_strict$outliers), nrow(lo_loose$outliers))
})

test_that("lookout old_version works", {
  lo <- lookout(X_clean, old_version = TRUE)
  expect_s3_class(lo, "lookoutliers")
})

test_that("lookout scale = FALSE works", {
  lo <- lookout(X_clean, scale = FALSE)
  expect_s3_class(lo, "lookoutliers")
})

test_that("lookout with pre-supplied gpd parameters works", {
  lo1 <- lookout(X_clean)
  lo2 <- lookout(X_clean, gpd = lo1$gpd)
  # Same GPD params should give same probabilities
  expect_equal(lo1$outlier_probability, lo2$outlier_probability)
})

test_that("lookout validates alpha range", {
  expect_error(lookout(X_clean, alpha = -0.1), "alpha")
  expect_error(lookout(X_clean, alpha = 1.1), "alpha")
})

test_that("lookout validates beta range", {
  expect_error(lookout(X_clean, beta = -0.1), "beta")
  expect_error(lookout(X_clean, beta = 1.5), "beta")
})

test_that("lookout validates gamma range", {
  expect_error(lookout(X_clean, gamma = -0.1), "gamma")
  expect_error(lookout(X_clean, gamma = 1.5), "gamma")
})

test_that("lookout works with univariate input", {
  X1d <- data.frame(x = c(rnorm(100), 10))
  lo <- lookout(X1d)
  expect_s3_class(lo, "lookoutliers")
})

# --- find_tda_bw() ---

test_that("find_tda_bw returns a positive scalar", {
  bw <- find_tda_bw(X_clean)
  expect_length(bw, 1L)
  expect_gt(bw, 0)
})

test_that("find_tda_bw warns about the deprecated fast argument", {
  expect_warning(bw <- find_tda_bw(X_clean, fast = TRUE), "deprecated")
  expect_equal(bw, find_tda_bw(X_clean))
  expect_warning(
    expect_equal(find_tda_bw(X_clean, fast = FALSE), find_tda_bw(X_clean)),
    "deprecated"
  )
})

test_that("find_tda_bw use_differences = TRUE works", {
  bw <- find_tda_bw(X_clean, use_differences = TRUE)
  expect_length(bw, 1L)
  expect_gt(bw, 0)
})

test_that("find_tda_bw validates gamma", {
  expect_error(find_tda_bw(X_clean, gamma = 0))
  expect_error(find_tda_bw(X_clean, gamma = 1.1))
})

# --- persisting_outliers() ---

test_that("persisting_outliers returns correct structure", {
  po <- persisting_outliers(X_outliers, num_steps = 5)
  expect_s3_class(po, "persistingoutliers")
  expect_equal(dim(po$out), c(NROW(X_outliers), 5L, 10L))
  expect_length(po$bw, 5L)
  expect_true(all(po$bw > 0))
  expect_true(!is.unsorted(po$bw))
})

test_that("persisting_outliers bandwidth grid is on the kernel support scale", {
  # Every element of `bw` is handed to lookout() as `bw`, so the whole grid
  # must be on the Epanechnikov support scale: the death radius multiplied by
  # sqrt(NCOL(X) + 4). The start of the grid used to be left unscaled.
  X <- mvscale(as.matrix(X_outliers))
  if (utils::packageVersion("mlpack") < "4.8.0") {
    death_radi <- mlpack::emst(X)$output[, 3]
  } else {
    death_radi <- mlpack::emst(X)[, 3]
  }
  kernel_scale <- sqrt(NCOL(X) + 4)

  po <- persisting_outliers(X, scale = FALSE, st_qq = 0.9, num_steps = 5)
  expect_equal(
    min(po$bw),
    unname(quantile(death_radi, probs = 0.9)) * kernel_scale
  )
  expect_equal(max(po$bw), max(death_radi) * kernel_scale)

  # st_qq = 0.9 is below lookout()'s gamma = 0.98, so the grid brackets the
  # bandwidth lookout() picks for itself.
  bw_lookout <- lookout(X, scale = FALSE)$bandwidth
  expect_lt(min(po$bw), bw_lookout)
  expect_gt(max(po$bw), bw_lookout)
})

test_that("persisting_outliers old_version works", {
  po <- persisting_outliers(X_outliers, num_steps = 5, old_version = TRUE)
  expect_s3_class(po, "persistingoutliers")
  expect_true(all(po$bw > 0))
})

# --- mvscale() ---

test_that("mvscale returns same dimensions as input", {
  z <- mvscale(X_clean)
  expect_equal(dim(z), dim(X_clean))
})

test_that("mvscale returns approximately unit variance", {
  z <- mvscale(X_clean)
  # After robust scaling columns have names z1, z2
  vars <- apply(as.matrix(z), 2, var)
  expect_true(all(vars > 0.1 & vars < 10))
})

test_that("mvscale works on a matrix", {
  z <- mvscale(as.matrix(X_clean), warning = FALSE)
  expect_true(is.matrix(z))
  expect_equal(dim(z), dim(X_clean))
})

test_that("mvscale cov = NULL skips rotation", {
  z <- mvscale(X_clean, cov = NULL, warning = FALSE)
  expect_equal(ncol(z), ncol(X_clean))
  # Column names unchanged when no rotation
  expect_equal(names(z), names(X_clean))
})

test_that("mvscale warns on non-numeric columns", {
  df <- cbind(X_clean, cat = letters[1:200])
  expect_warning(mvscale(df), "non-numeric")
})

test_that("mvscale errors on non-numeric vector", {
  expect_error(mvscale(letters[1:10]))
})

test_that("mvscale returns a vector with terms attached for vector input", {
  z <- mvscale(X_clean$x)
  expect_type(z, "double")
  expect_length(z, NROW(X_clean))
  expect_equal(median(z), 0, tolerance = 1e-6)
  expect_equal(attr(z, "center"), median(X_clean$x))
  expect_equal(attr(z, "scale"), robustbase::s_Qn(X_clean$x))
  expect_equal(attr(z, "scale_inverse"), 1 / attr(z, "scale"))
})

test_that("mvscale returns the covariance and its inverse as terms", {
  z <- mvscale(X_clean)
  S <- attr(z, "scale")
  Sinv <- attr(z, "scale_inverse")
  expect_equal(dim(S), c(2L, 2L))
  expect_equal(unname(S %*% Sinv), diag(2), tolerance = 1e-6)
  expect_equal(attr(z, "center"), apply(X_clean, 2, median))
})

test_that("mvscale cov = NULL returns the per-column scales as terms", {
  z <- mvscale(X_clean, cov = NULL)
  expect_equal(attr(z, "scale"), apply(X_clean, 2, robustbase::s_Qn))
  expect_equal(attr(z, "scale_inverse"), 1 / attr(z, "scale"))
})

test_that("mvscale renames matrix columns only when rotation is applied", {
  z <- mvscale(as.matrix(X_clean))
  expect_equal(colnames(z), c("z1", "z2"))
  expect_null(names(z))
  z <- mvscale(as.matrix(X_clean), cov = NULL)
  expect_equal(colnames(z), names(X_clean))
})

test_that("mvscale does not rename a single column", {
  expect_equal(names(mvscale(X_clean["x"])), "x")
})

test_that("mvscale passes alpha and extra arguments to cov()", {
  z <- mvscale(X_clean, alpha = 0.5)
  expect_equal(dim(z), dim(X_clean))
  expect_false(isTRUE(all.equal(attr(z, "scale"), attr(mvscale(X_clean), "scale"))))
  z <- mvscale(X_clean, cov = stats::cov, use = "complete.obs")
  expect_equal(attr(z, "scale"), stats::cov(as.matrix(X_clean)))
})

test_that("mvscale omits missing values when estimating terms", {
  X_na <- X_clean
  X_na[3, "x"] <- NA
  expect_equal(
    attr(mvscale(X_na, cov = NULL), "center"),
    apply(X_na, 2, median, na.rm = TRUE)
  )
  # The rotation makes the whole row missing
  z <- mvscale(X_na)
  expect_true(all(is.na(z[3, ])))
  expect_false(anyNA(z[-3, ]))
})

test_that("mvscale errors on infinite values", {
  X_inf <- X_clean
  X_inf[1, "y"] <- Inf
  expect_error(mvscale(X_inf), "infinite")
  expect_error(mvscale(c(X_clean$x, -Inf)), "infinite")
})

test_that("mvscale errors on non-numeric matrices and unsupported objects", {
  expect_error(mvscale(matrix(letters[1:10], ncol = 2)), "numeric")
  expect_error(mvscale(array(1:24, c(2, 3, 4))), "vector, matrix or data frame")
})

test_that("mvscale errors when cov() returns no covariance matrix", {
  expect_error(mvscale(X_clean, cov = function(x, ...) "rubbish"))
  expect_error(mvscale(X_clean, cov = function(x, ...) list(a = 1)), "can't find")
})

test_that("lookout uses the fast approximation by default only when n > 10000", {
  fast_default <- formals(lookout)$fast
  expect_false(eval(fast_default, list(X = matrix(0, 10000, 2))))
  expect_true(eval(fast_default, list(X = matrix(0, 10001, 2))))
})

test_that("lookout explains why the GPD cannot be fitted when too many points are isolated", {
  # 700 tightly clustered points and 300 scattered ones. With bw = 0.5 almost
  # all of the scattered points have no other point inside the kernel support,
  # so well over 20% of surprisals tie at their maximum, which is more than
  # 100(1 - beta)% for both beta = 0.9 and beta = 0.8. The fast path keeps
  # every point's nearest neighbours, so it cannot make a point isolated and
  # fails in exactly the same way as the exact path.
  set.seed(5)
  X <- matrix(c(rnorm(700, sd = 0.001), runif(300, 0, 2000)), ncol = 1)
  expect_error(
    lookout(X, bw = 0.5, scale = FALSE, fast = FALSE),
    "more than 10% of observations have no other observation inside the kernel support",
    fixed = TRUE
  )
  expect_error(
    lookout(X, bw = 0.5, scale = FALSE, fast = TRUE),
    "more than 10% of observations have no other observation inside the kernel support",
    fixed = TRUE
  )
  expect_error(lookout(X, bw = 0.5, scale = FALSE, beta = 0.8), "more than 20%")
  expect_error(lookout(X, bw = 0.5, scale = FALSE), "larger `gamma`", fixed = TRUE)
  # A wider kernel support removes the problem
  expect_s3_class(lookout(X, bw = 30, scale = FALSE), "lookoutliers")
})
