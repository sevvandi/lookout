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
#>    Outliers Probability
#> 1        27           0
#> 2        49           0
#> 5        57           0
#> 6        65           0
#> 7       114           0
#> 8       119           0
#> 9       132           0
#> 10      140           0
#> 11      153           0
#> 12      177           0
#> 14      193           0
#> 
```
