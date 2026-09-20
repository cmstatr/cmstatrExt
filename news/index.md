# Changelog

## Version 0.5.0

CRAN release: 2026-09-20

- Added function for computing “one-sample” equivalency factors. This
  function should produce the same values as
  [`cmstatr::k_equiv()`](https://www.cmstatr.net/reference/k_equiv.html),
  but is written in C++ and is hence faster.
- Renamed the function
  [`p_equiv()`](https://cmstatrExt.cmstatr.net/reference/p_equiv.md) to
  [`p_equiv_one_sample()`](https://cmstatrExt.cmstatr.net/reference/p_equiv_one_sample.md).
  The function
  [`p_equiv()`](https://cmstatrExt.cmstatr.net/reference/p_equiv.md)
  will still work but will produce a warning and may be removed in a
  future version of this package.
- Reduced the tolerance for root finding from `DBL_EPSILON^0.25` (about
  1.2e-4) `DBL_EPSILON^0.5` (about 1.5e-8). This may affect the results
  of `k_equiv_two_sample`, `p_equiv_two_sample` and `p_equiv_one_sample`
  and should produce more accurate results for each of these functions.
  You can expect to see the largest differences with very low p-values
  and alpha values.
- In the function
  [`power_sim_dual()`](https://cmstatrExt.cmstatr.net/reference/power_sim_dual.md),
  changed the internal data types used to count the number of
  acceptance/equivalency failures from `int` to `unsigned long long` to
  prevent overflow when `replicates` is greater than
  46340. 
- This package now checks for conflicting function names and will
  produce a message when attached, similar to the way that the
  `tidyverse` package produces a conflict message when attached.
- In the accompanying website \[<https://cmstatrExt.cmstatr.net>\],
  added a new calculator for equivalency thresholds.

## Version 0.4.1

CRAN release: 2026-03-18

- Minor update to call Rcpp::stop in lieu of Rf_error due to an upcoming
  change to the Rcpp package. No changed in functionality.

## Version 0.4.0

CRAN release: 2024-05-07

- First release on CRAN
- Minor documentation improvements

## Version 0.3.0

- Added functions for creating average stress-strain curves
  (`average_curve_lm` and `average_curve_optim`)
- Added example stress-strain data (`pa12_tension` and `fff_shear`)
- Added vignette with examples of fitting average curves to
  stress-strain data

## Version 0.2.1

- Update to p-value vignette

## Version 0.2.0

- Created p-value calculator for website
- Added `iso_equiv_two_sample` function
- Created p-value vignette
- Added more unit and integration tests

## Version 0.1.0

- Documentation improvements
- Added equivalency/acceptance factor calculator to website

## Version 0.1.0

- Added power simulation function
- Added documentation website
- Changed c++ unit test framework to testthat (which uses catch2)
- Improved speed of numerical integration

## Version 0.0.1

- First public release
