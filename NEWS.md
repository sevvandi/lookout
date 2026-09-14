# lookout (development version)

- Fixed the Epanechnikov kernel scaling for multivariate data, which affected the leave-one-out density estimates when `NCOL(X) > 1`.
- `mvscale()` now matches `weird::mvscale()`: the default covariance estimate is a robust MCD estimate (with new `alpha` argument) rather than OGK, extra arguments are passed to `cov()`, missing values are omitted when estimating the center, scale and covariance, infinite values throw an error, and the center, scale and inverse scale are returned as attributes.
  This changes the scaled data used by `lookout()` and `persisting_outliers()` when `scale = TRUE`.
- The `fast` argument of find_tda_bw() is deprecated and ignored.
- Bug fixes and documentation improvements

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
