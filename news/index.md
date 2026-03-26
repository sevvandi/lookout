# Changelog

## lookout (development version)

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
