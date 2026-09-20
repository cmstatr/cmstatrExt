# Calculate the factors for a one-sample acceptance test

Calculates the factors k1 and k2, which are used for setting acceptance
values for lot acceptance. These factors consider the size of acceptance
sample (`m`).

## Usage

``` r
k_equiv_one_sample(alpha, m)
```

## Arguments

- alpha:

  the nominal significance of the test

- m:

  the size of the acceptance sample

## Value

A vector of length 2 with the contents `c(k1, k2)`

## Details

This function is equivalent to
[`cmstatr::k_equiv()`](https://www.cmstatr.net/reference/k_equiv.html),
but is implemented in C++ instead of in R, and is hence slightly faster.

## References

Vangel, M. (2002). Lot Acceptance and Compliance Testing Using the
Sample Mean and an Extremum. Technometrics.
https://doi.org/10.1198/004017002188618428
