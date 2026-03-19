# Augment a `data.frame` with the results from `average_curve_lm`

Augment a `data.frame` with the results from `average_curve_lm`

## Usage

``` r
# S3 method for class 'average_curve_lm'
augment(x, newdata = NULL, extrapolate = FALSE, ...)
```

## Arguments

- x:

  an `average_curve_lm` object

- newdata:

  (optional) a new `data.frame` to which to augment the object

- extrapolate:

  whether to show the curve fit on all data or only the data within the
  original fitted range. Default: FALSE

- ...:

  ignored

## Value

a `data.frame` with new columns `.fit`, `.extrapolate` and `.residual`

## See also

[`average_curve_lm()`](https://cmstatrExt.cmstatr.net/reference/average_curve_lm.md)

## Examples

``` r
curve_fit <- average_curve_lm(
  pa12_tension,
  Coupon,
  Stress ~ I(Strain) + I(Strain^2) + I(Strain^3) + 0,
  n_bins = 100
)
augment(curve_fit)
#> # A tibble: 3,105 × 6
#>    Coupon     Strain Stress  .fit .extrapolate .residual
#>    <chr>       <dbl>  <dbl> <dbl> <lgl>            <dbl>
#>  1 Coupon 4 0        0.0561 0     FALSE           0.0561
#>  2 Coupon 4 0.000200 0.247  0.235 FALSE           0.0122
#>  3 Coupon 4 0.000400 0.569  0.468 FALSE           0.100 
#>  4 Coupon 4 0.000601 0.440  0.702 FALSE          -0.262 
#>  5 Coupon 4 0.000801 0.778  0.934 FALSE          -0.156 
#>  6 Coupon 4 0.00100  0.854  1.17  FALSE          -0.312 
#>  7 Coupon 4 0.00120  0.955  1.40  FALSE          -0.442 
#>  8 Coupon 4 0.00140  1.40   1.63  FALSE          -0.230 
#>  9 Coupon 4 0.00160  1.54   1.86  FALSE          -0.320 
#> 10 Coupon 4 0.00180  1.62   2.09  FALSE          -0.463 
#> # ℹ 3,095 more rows
## # A tibble: 3,105 × 6
##    Coupon     Strain  Stress  .fit .extrapolate .residual
##    <chr>       <dbl>   <dbl> <dbl> <lgl>            <dbl>
##  1 Coupon 4 0        -0.353  0     FALSE          -0.353
##  2 Coupon 4 0.000200 -0.0604 0.235 FALSE          -0.295
##  3 Coupon 4 0.000400  0.283  0.469 FALSE          -0.185
##  4 Coupon 4 0.000601  0.475  0.702 FALSE          -0.228
##  5 Coupon 4 0.000801  0.737  0.935 FALSE          -0.198
##  6 Coupon 4 0.00100   0.803  1.17  FALSE          -0.364
##  7 Coupon 4 0.00120   1.25   1.40  FALSE          -0.151
##  8 Coupon 4 0.00140   1.32   1.63  FALSE          -0.305
##  9 Coupon 4 0.00160   1.53   1.86  FALSE          -0.325
## 10 Coupon 4 0.00180   2.01   2.09  FALSE          -0.0735
## # i 3,095 more row
## # i Use `print(n = ...)` to see more rows
```
