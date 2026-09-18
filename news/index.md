# Changelog

## lookout (development version)

- The exact density estimates (`fast = FALSE`) now use a fixed-radius
  neighbour search via
  [`dbscan::frNN()`](https://rdrr.io/pkg/dbscan/man/frNN.html), so
  memory no longer grows with `NROW(X)^2`. Results are unchanged up to
  rounding.
- Nearest neighbour searches now use dbscan rather than RANN, which is
  no longer imported.
- The default is now `fast = NROW(X) > 10000`.
- [`lookout()`](https://sevvandi.github.io/lookout/reference/lookout.md)
  gives a more helpful error when the GPD cannot be fitted because too
  many observations have no neighbours within the bandwidth.
- [`mvscale()`](https://sevvandi.github.io/lookout/reference/mvscale.md)
  now matches `weird::mvscale()`, with a robust MCD covariance estimate
  by default. This changes the scaled data used by
  [`lookout()`](https://sevvandi.github.io/lookout/reference/lookout.md)
  and
  [`persisting_outliers()`](https://sevvandi.github.io/lookout/reference/persisting_outliers.md)
  when `scale = TRUE`.
- Fixed the Epanechnikov kernel scaling for multivariate data.
- The `fast` argument of
  [`find_tda_bw()`](https://sevvandi.github.io/lookout/reference/find_tda_bw.md)
  is deprecated and ignored.
- Now requires R \>= 4.1.0.
- Bug fixes and documentation improvements.

## lookout 2.0.2

CRAN release: 2026-07-21

- Compatibility with mlpack 4.8.0, which changed the value returned by
  [`mlpack::emst()`](https://rdrr.io/pkg/mlpack/man/emst.html).

## lookout 2.0.1

CRAN release: 2026-03-26

- find_tda_bw() now much faster, especially for large data sets, by
  using MST instead of Rips complex to compute persistent homology.

## lookout 2.0.0

CRAN release: 2026-01-19

- Updated lookout algorithm as per Hyndman, Kandanaarachchi and Turner
  (2026).
- Added mvscale() to do robust multivariate scaling.
- Exported find_tda_bw().

## lookout 0.1.4

CRAN release: 2022-10-13

- First CRAN version, based on Kandanaarachchi and Hyndman (JCGS, 2022).
