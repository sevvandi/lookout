# Compute robust multivariate scaled data

A multivariate version of
[`base::scale()`](https://rdrr.io/r/base/scale.html), that takes account
of the covariance matrix of the data. By default, robust estimates are
used: the centers are removed using medians, the scale function for
univariate data is `s_Qn`, and the covariance matrix for multivariate
data is estimated using a robust MCD estimate. The data are scaled using
the Cholesky decomposition of the inverse (co)variance. Then the scaled
data are returned. Details of the methods are provided by Hyndman
(2026).

## Usage

``` r
mvscale(
  object,
  center = stats::median,
  scale = robustbase::s_Qn,
  cov = robustbase::covMcd,
  alpha = 0.9,
  warning = TRUE,
  ...
)
```

## Arguments

- object:

  A vector, matrix, or data frame containing some numerical data.

- center:

  A function to compute the center of each numerical variable. Set to
  NULL if no centering is required.

- scale:

  A function to scale each numerical variable. When
  `cov = robustbase::covOGK()`, `scale` is passed as the `sigmamu`
  argument. When `cov = robustbase::covMcd()`, `scale` is passed as the
  `scalefn` argument.

- cov:

  A function to compute the covariance matrix. Set to NULL if no
  rotation required. [`cov()`](https://rdrr.io/r/stats/cor.html) must
  either return the matrix directly, or a list containing a matrix named
  `cov`.

- alpha:

  When `cov = robustbase::covMcd()`, `alpha` controls the size of the
  subsets over which the determinant is minimized. Otherwise it is
  ignored. Set to 0.9 by default.

- warning:

  Should a warning be issued if non-numeric columns are ignored?

- ...:

  Other arguments are passed to
  [`cov()`](https://rdrr.io/r/stats/cor.html).

## Value

A vector, matrix or data frame of the same size and class as `object`,
but with numerical variables replaced by scaled versions (renamed if
they have been rotated).

## Details

Optionally, the centering and scaling can be done for each variable
separately, by setting `cov = NULL`, so there is no rotation of the
data, Also optionally, non-robust methods can be used by specifying
`center = mean`, `scale = stats::sd()`, and `cov = stats::cov()`. Any
non-numeric columns are retained with a warning. Missing values are
removed before the centers, scale and cov are estimated.

## References

Hyndman, R J (2026) "That's weird: Anomaly detection using R", Section
3.7.

## See also

[`base::scale()`](https://rdrr.io/r/base/scale.html),
[`stats::sd()`](https://rdrr.io/r/stats/sd.html),
[`stats::cov()`](https://rdrr.io/r/stats/cor.html),
[`robustbase::covMcd()`](https://rdrr.io/pkg/robustbase/man/covMcd.html),
[`robustbase::covOGK()`](https://rdrr.io/pkg/robustbase/man/covOGK.html),
[`robustbase::s_Qn()`](https://rdrr.io/pkg/robustbase/man/Qn.html)

## Author

Rob J Hyndman

## Examples

``` r
# Univariate z-scores
z <- mvscale(faithful$eruptions, center = mean, scale = sd)
# Non-robust scaling with no rotation
z <- mvscale(faithful, center = mean, scale = sd, cov = NULL)
# Non-robust scaling with rotation
z <- mvscale(faithful, center = mean, scale = sd, cov = stats::cov)
# Robust scaling and rotation
z <- mvscale(faithful)
```
