# Identifies outliers using the algorithm lookout.

This function identifies outliers using the algorithm lookout, an
outlier detection method that uses leave-one-out kernel density
estimates and generalized Pareto distributions to find outliers.

## Usage

``` r
lookout(
  X,
  alpha = 0.01,
  beta = 0.9,
  gamma = 0.98,
  bw = NULL,
  gpd = NULL,
  scale = TRUE,
  fast = NROW(X) > 1e+05,
  old_version = FALSE
)
```

## Arguments

- X:

  The numerical input data in a data.frame, matrix or tibble format.

- alpha:

  The level of significance. Default is `0.01`. So there is a 1/100
  chance of any point being falsely classified as an outlier.

- beta:

  The quantile threshold used in the GPD estimation. Default is `0.90`.
  To ensure there is enough data available, values greater than 0.90 are
  set to 0.90.

- gamma:

  Parameter for bandwidth calculation giving the quantile of the Rips
  death radii to use for the bandwidth. Default is `0.98`. Ignored under
  the old version; where the lower limit of the maximum Rips death radii
  difference is used. Also ignored if `bw` is provided.

- bw:

  Bandwidth parameter. If `NULL` (default), the bandwidth is found using
  Persistent Homology.

- gpd:

  Generalized Pareto distribution parameters. If `NULL` (the default),
  these are estimated from the data.

- scale:

  If `TRUE`, the data is standardized. Using the old version, unit
  scaling is applied so that each column is in the range `[0,1]`. Under
  the new version, robust rotation and scaling is used so that the
  columns are approximately uncorrelated with unit variance. Default is
  `TRUE`.

- fast:

  If `TRUE`, each density estimate uses only the `k` nearest neighbours
  of the point, where `k` is between 100 and 500 depending on `NROW(X)`,
  rather than all observations. This is an approximation: it is exact
  only when fewer than `k` observations lie within `bw` of every point,
  which is often false in more than two or three dimensions. Default is
  `TRUE` when `NROW(X) > 100000`. The bandwidth calculation always uses
  all of the data.

- old_version:

  Logical indicator of which version of the algorithm to use. Default is
  FALSE, meaning the newer version is used.

## Value

A list with the following components:

- `data`:

  The input data `X`, before any scaling.

- `outliers`:

  The set of outliers.

- `outlier_probability`:

  The GPD probability of the data.

- `outlier_scores`:

  The outlier scores of the data.

- `bandwidth`:

  The bandwdith selected using persistent homology.

- `kde`:

  The kernel density estimate values.

- `lookde`:

  The leave-one-out kde values.

- `gpd`:

  The fitted GPD parameters.

- `call`:

  The matched call.

## References

Kandanaarachchi, S, and Hyndman, RJ (2022) Leave-one-out kernel density
estimates for outlier detection, *J Computational & Graphical
Statistics*, **31**(2), 586-599.
<https://robjhyndman.com/publications/lookout/>.

Hyndman, RJ, Kandanaarachchi, S, and Turner, K (2026) When lookout meets
crackle: Anomaly detection using kernel density estimation, in
preparation. <https://robjhyndman.com/publications/lookout2.html>

## Examples

``` r
X <- rbind(
  data.frame(
    x = rnorm(500),
    y = rnorm(500)
  ),
  data.frame(
    x = rnorm(5, mean = 10, sd = 0.2),
    y = rnorm(5, mean = 10, sd = 0.2)
  )
)
lo <- lookout(X)
lo
#> Leave-out-out KDE outliers using lookout algorithm
#> 
#> Call: lookout(X = X)
#> 
#>    Outliers  Probability
#> 1       105 0.0026920212
#> 2       213 0.0006154015
#> 3       220 0.0051289400
#> 4       298 0.0000000000
#> 5       310 0.0017968741
#> 6       472 0.0007241131
#> 7       501 0.0094663929
#> 8       502 0.0094780204
#> 9       503 0.0078094157
#> 10      504 0.0097237152
#> 11      505 0.0097534626
#> 
autoplot(lo)
```
