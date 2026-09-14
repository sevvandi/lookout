# Identifies bandwidth for outlier detection.

This function identifies the bandwidth that is used in the kernel
density estimate computation. The function uses topological data
analysis (TDA) to find the badnwidth.

## Usage

``` r
find_tda_bw(X, fast = NULL, gamma = 0.98, use_differences = FALSE)
```

## Arguments

- X:

  The numerical input data in a data.frame, matrix or tibble format.

- fast:

  Deprecated and ignored. The bandwidth is always computed using all of
  `X`.

- gamma:

  Parameter for bandwidth calculation giving the quantile of the Rips
  death radii to use for the bandwidth. Default is `0.98`. Ignored under
  the old version; where the lower limit of the maximum Rips death radii
  difference is used. Also ignored if `bw` is provided.

- use_differences:

  If TRUE, the bandwidth is set to the lower point of the maximum Rips
  death radii differences. If FALSE, the gamma quantile of the Rips
  death radii is used. Default is FALSE.

## Value

The bandwidth

## Details

The value returned is the raw quantile of the minimum spanning tree edge
lengths.
[`lookout`](https://sevvandi.github.io/lookout/reference/lookout.md)
multiplies it by `sqrt(m + 4)`, where `m = NCOL(X)`, to obtain the
support radius of the Epanechnikov kernel, so that the kernel has
marginal standard deviation equal to this quantile.

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
find_tda_bw(X)
#> [1] 0.4605964
```
