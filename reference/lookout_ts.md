# Identifies outliers in univariate time series using the algorithm lookout.

This is the time series implementation of lookout which identifies
outliers in the double differenced time series.

## Usage

``` r
lookout_ts(x, scale = FALSE, ...)
```

## Arguments

- x:

  The input univariate time series.

- scale:

  If `TRUE`, the data is standardized. Using the old version, unit
  scaling is applied so that each column is in the range `[0,1]`. Under
  the new version, robust rotation and scaling is used so that the
  columns are approximately uncorrelated with unit variance. Default is
  `TRUE`.

- ...:

  Other arguments are passed to
  [`lookout`](https://sevvandi.github.io/lookout/reference/lookout.md).

## Value

A lookout object.

## See also

[`lookout`](https://sevvandi.github.io/lookout/reference/lookout.md)

## Examples

``` r
set.seed(1)
x <- arima.sim(list(order = c(1, 1, 0), ar = 0.8), n = 200)
x[50] <- x[50] + 10
plot(x)

lo <- lookout_ts(x)
lo
#> Leave-out-out KDE outliers using lookout algorithm
#> 
#> Call: lookout(X = u, scale = scale)
#> 
#>   Outliers Probability
#> 2       50 0.000000000
#> 4      177 0.006841116
#> 
```
