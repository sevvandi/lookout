# Computes the bandwidth for lookout from the minimum spanning tree.

This function computes the bandwidth used by
[`lookout()`](https://sevvandi.github.io/lookout/reference/lookout.md)
for its kernel density estimates, from the edge lengths of the Euclidean
minimum spanning tree of the data.

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

  The quantile of the minimum spanning tree edge lengths to return.
  Default is `0.98`. Ignored when `use_differences = TRUE`.

- use_differences:

  If `TRUE`, the value returned is the lower end of the largest gap
  between consecutive minimum spanning tree edge lengths among those at
  or above their median. If `FALSE` (the default), the `gamma` quantile
  of the edge lengths is returned.

## Value

A single number: the bandwidth on the scale of `X`, before
multiplication by `sqrt(NCOL(X) + 4)`.

## Details

The Euclidean minimum spanning tree of the rows of `X` is computed with
[`emst()`](https://rdrr.io/pkg/mlpack/man/emst.html), giving
`NROW(X) - 1` edge lengths in increasing order. Each edge length is the
distance at which two clusters merge under single-linkage clustering, or
equivalently the death radius of a connected component in the Rips
filtration of the data, which is why the function is named after
topological data analysis.

By default (`use_differences = FALSE`), the value returned is the
`gamma` quantile of the edge lengths, computed with
[`quantile()`](https://rdrr.io/r/stats/quantile.html) using `type = 8`.
With `use_differences = TRUE`, only the edge lengths at or above their
median are kept, and the value returned is the lower of the two
consecutive edge lengths with the largest gap between them. This is the
rule of Kandanaarachchi and Hyndman (2022), used by
[`lookout()`](https://sevvandi.github.io/lookout/reference/lookout.md)
when `old_version = TRUE`.

The value returned is on the scale of `X` and is not itself the support
radius of the kernel.
[`lookout()`](https://sevvandi.github.io/lookout/reference/lookout.md)
multiplies it by `sqrt(m + 4)`, where `m = NCOL(X)`, to obtain the
support radius of the Epanechnikov kernel. A spherically symmetric
Epanechnikov kernel with support radius `h` in `m` dimensions has
standard deviation `h / sqrt(m + 4)` in each coordinate, so the kernel
used by
[`lookout()`](https://sevvandi.github.io/lookout/reference/lookout.md)
has per-coordinate standard deviation equal to the value returned here.
Note that
[`lookout()`](https://sevvandi.github.io/lookout/reference/lookout.md)
applies this function to the scaled data when `scale = TRUE`.

The minimum spanning tree is always computed on all of `X`.

## References

Kandanaarachchi, S, and Hyndman, RJ (2022) Leave-one-out kernel density
estimates for outlier detection, *J Computational & Graphical
Statistics*, **31**(2), 586-599.
<https://robjhyndman.com/publications/lookout/>.

Hyndman, RJ, Kandanaarachchi, S, and Turner, K (2026) Lookout 2: Anomaly
detection via leave-one-out kernel density estimation, arXiv:2603.22636.
<https://robjhyndman.com/publications/lookout2.html>

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
