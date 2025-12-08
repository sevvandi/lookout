# Identifies bandwidth for outlier detection.

This function identifies the bandwidth that is used in the kernel
density estimate computation. The function uses topological data
analysis (TDA) to find the badnwidth.

## Usage

``` r
find_tda_bw(X, fast = TRUE, gamma = 0.97, use_differences = FALSE)
```

## Arguments

- X:

  The numerical input data in a data.frame, matrix or tibble format.

- fast:

  If `TRUE` (default), makes the computation faster by sub-setting the
  data for the bandwidth calculation.

- gamma:

  Parameter for bandwidth calculation giving the quantile of the Rips
  death radii to use for the bandwidth. Default is `0.97`. Ignored under
  the old version; where the lower limit of the maximum Rips death radii
  difference is used. Also ignored if `bw` is provided.

- use_differences:

  If TRUE, the bandwidth is set to the lower point of the maximum Rips
  death radii differences. If FALSE, the gamma quantile of the Rips
  death radii is used. Default is FALSE.

## Value

The bandwidth

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
find_tda_bw(X, fast = TRUE)
#> [1] 2.222265
```
