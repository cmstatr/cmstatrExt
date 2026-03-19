# Generate an average curve using `lm`

The user must decide on a single dependent variable (`Y`) and a single
independent variable (`X`). The user will specify a `formula` with the
relationship between the dependent and independent variables. For a
`data.frame` containing stress-strain (or load-deflection) data for more
than one coupon, the maximum value of `X` for each coupon is found and
the smallest maximum value determines the range over which the curve fit
is performed: the range is from zero to this value. Only positive values
of `X` are considered. For each coupon individually, the data is divided
into a user-specified number of bins and averaged within each bin. The
resulting binned/averaged data is then passed to
[`stats::lm()`](https://rdrr.io/r/stats/lm.html) to perform the curve
fitting.

## Usage

``` r
average_curve_lm(data, coupon_var, model, n_bins = 100)
```

## Arguments

- data:

  a `data.frame`

- coupon_var:

  the variable for coupon identification

- model:

  a `formula` for the curve to fit

- n_bins:

  the number of bins to average the data inside into before fitting

## Value

an object of class `average_curve_lm` with the following content:

- `data` the original data provided to the function

- `binned_data` the data after the binning/averaging operation

- `fit_lm` the results of the call to `lm`

- `n_bins` the number of bins specified by the user

- `max_x` the upper end of the range used for fitting

- `y_var` the independent (`Y`) variable

- `x_var` the dependent (`X`) variable

## Details

When specifying the formula (argument `model`), there are two things to
keep in mind. First, based on physical behavior, it is normally
desirable to set the intercept to zero (e.g. so that there is 0 stress
at 0 strain). To do this, include a term `+0` in the formula. Second,
when specifying a term for a power of the `X` variable (for example,
\$X^2\$), this needs to be wrapped inside the "as-is" operator
[`I()`](https://rdrr.io/r/base/AsIs.html), otherwise, `R` will treat it
as an interaction term, rather than an exponent. In other words, if you
want to include a quadratic term, you need to write `I(X^2)` (replacing
`X` with the appropriate variable from your `data.frame`).

## See also

[`~`](https://rdrr.io/r/base/tilde.html),
[`I()`](https://rdrr.io/r/base/AsIs.html),
[`lm()`](https://rdrr.io/r/stats/lm.html),
[`average_curve_optim()`](https://cmstatrExt.cmstatr.net/reference/average_curve_optim.md),
[`print.average_curve_lm()`](https://cmstatrExt.cmstatr.net/reference/print.average_curve_lm.md),
[`summary.average_curve_lm()`](https://cmstatrExt.cmstatr.net/reference/summary.average_curve_lm.md),
[`augment.average_curve_lm()`](https://cmstatrExt.cmstatr.net/reference/augment.average_curve_lm.md)

## Examples

``` r
# using the `pa12_tension` dataset and fitting a cubic polynomial with
# zero intercept:
curve_fit <- average_curve_lm(
  pa12_tension,
  Coupon,
  Stress ~ I(Strain) + I(Strain^2) + I(Strain^3) + 0,
  n_bins = 100
)
print(curve_fit)
#> 
#> Range: ` Strain ` in  [ 0,  0.1409409 ]
#> 
#> Call:
#> average_curve_lm(data = pa12_tension, coupon_var = Coupon, model = Stress ~ 
#>     I(Strain) + I(Strain^2) + I(Strain^3) + 0, n_bins = 100)
#> 
#> Coefficients:
#>   I(Strain)  I(Strain^2)  I(Strain^3)  
#>        1173        -8762        20481  
#> 
## Range: ` Strain ` in  [ 0,  0.1409409 ]
##
## Call:
##   average_curve_lm(data = pa12_tension, coupon_var = Coupon,
##                    model = Stress ~ I(Strain) + I(Strain^2) + I(Strain^3)
##                    + 0, n_bins = 100)
##
## Coefficients:
##    I(Strain)   I(Strain^2)   I(Strain^3)
##        1174         -8783         20586
```
