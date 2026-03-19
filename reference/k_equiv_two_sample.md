# Calculate the factors for a two-sample acceptance test

Calculates the factors k1 and k2, which are used for setting acceptance
values for lot acceptance. These factors consider both the size of the
qualification sample (`n`) and the size of acceptance sample (`m`). This
test is detailed in a forthcoming paper.

## Usage

``` r
k_equiv_two_sample(alpha, n, m)
```

## Arguments

- alpha:

  the desired probability of Type 1 error

- n:

  the size of the qualification sample

- m:

  the size of the acceptance sample

## Value

A vector of length 2 with the contents `c(k1, k2)`

## References

Kloppenborg, S. (2023). Lot acceptance testing using sample mean and
extremum with finite qualification samples. Journal of Quality
Technology, https://doi.org/10.1080/00224065.2022.2147884
