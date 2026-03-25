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

  If `TRUE` (default), makes the computation faster by sub-setting the
  data for the bandwidth calculation.

- old_version:

  Logical indicator of which version of the algorithm to use. Default is
  FALSE, meaning the newer version is used.

## Value

A list with the following components:

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
#>   Outliers  Probability
#> 1      105 1.710857e-03
#> 2      213 2.880705e-04
#> 3      220 4.790137e-03
#> 4      298 0.000000e+00
#> 5      310 1.088849e-06
#> 6      329 9.808537e-03
#> 7      351 9.915623e-03
#> 8      472 1.820545e-04
#> 9      503 9.930251e-03
#> 
autoplot(lo)
```
