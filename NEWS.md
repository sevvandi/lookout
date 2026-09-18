# lookout (development version)

- The exact density estimates (`fast = FALSE`) now use a fixed-radius neighbour search via `dbscan::frNN()`, so memory no longer grows with `NROW(X)^2`. Results are unchanged up to rounding.
- Nearest neighbour searches now use dbscan rather than RANN, which is no longer imported.
- The default is now `fast = NROW(X) > 10000`.
- `lookout()` gives a more helpful error when the GPD cannot be fitted because too many observations have no neighbours within the bandwidth.
- `mvscale()` now matches `weird::mvscale()`, with a robust MCD covariance estimate by default. This changes the scaled data used by `lookout()` and `persisting_outliers()` when `scale = TRUE`.
- Fixed the Epanechnikov kernel scaling for multivariate data.
- The `fast` argument of `find_tda_bw()` is deprecated and ignored.
- Now requires R >= 4.1.0.
- Bug fixes and documentation improvements.

# lookout 2.0.2

- Compatibility with mlpack 4.8.0, which changed the value returned by `mlpack::emst()`.

# lookout 2.0.1

- find_tda_bw() now much faster, especially for large data sets, by using MST instead of Rips complex to compute persistent homology.

# lookout 2.0.0

- Updated lookout algorithm as per Hyndman, Kandanaarachchi and Turner (2026).
- Added mvscale() to do robust multivariate scaling.
- Exported find_tda_bw().

# lookout 0.1.4

- First CRAN version, based on Kandanaarachchi and Hyndman (JCGS, 2022).
