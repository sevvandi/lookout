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
  expect_error(lookout(X_clean, alpha = -0.1))
  expect_error(lookout(X_clean, alpha = 1.1))
})

test_that("lookout validates beta range", {
  expect_error(lookout(X_clean, beta = -0.1))
  expect_error(lookout(X_clean, beta = 1.5))
})

test_that("lookout validates gamma range", {
  expect_error(lookout(X_clean, gamma = -0.1))
  expect_error(lookout(X_clean, gamma = 1.5))
})

test_that("lookout works with univariate input", {
  X1d <- data.frame(x = c(rnorm(100), 10))
  lo <- lookout(X1d)
  expect_s3_class(lo, "lookoutliers")
})

# --- find_tda_bw() ---

test_that("find_tda_bw returns a positive scalar", {
  bw <- find_tda_bw(X_clean, fast = TRUE)
  expect_length(bw, 1L)
  expect_gt(bw, 0)
})

test_that("find_tda_bw fast and slow give similar results", {
  bw_fast <- find_tda_bw(X_clean, fast = TRUE)
  bw_slow <- find_tda_bw(X_clean, fast = FALSE)
  # Should be within an order of magnitude
  expect_true(abs(log(bw_fast) - log(bw_slow)) < log(10))
})

test_that("find_tda_bw use_differences = TRUE works", {
  bw <- find_tda_bw(X_clean, fast = TRUE, use_differences = TRUE)
  expect_length(bw, 1L)
  expect_gt(bw, 0)
})

test_that("find_tda_bw validates gamma", {
  expect_error(find_tda_bw(X_clean, gamma = 0))
  expect_error(find_tda_bw(X_clean, gamma = 1.1))
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
