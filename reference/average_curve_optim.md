# Generate an average curve using `optim`

The user must decide on a single dependent variable (`Y`) and a single
independent variable (`X`). The user will specify a function defining
the relationship between the dependent and independent variables. For a
`data.frame` containing stress-strain (or load-deflection) data for more
than one coupon, the maximum value of `X` for each coupon is found and
the smallest maximum value determines the range over which the curve fit
is performed: the range is from zero to this value. Only positive values
of `X` are considered. For each coupon individually, the data is divided
into a user-specified number of bins and averaged within each bin. The
resulting binned/averaged data is then used for curve fitting. The mean
squared error between the observed value of `Y` and the result of the
user-specified function evaluated at each `X` is minimized by varying
the parameters `par`.

## Usage

``` r
average_curve_optim(
  data,
  coupon_var,
  x_var,
  y_var,
  fn,
  par,
  n_bins = 100,
  method = "L-BFGS-B",
  ...
)
```

## Arguments

- data:

  a `data.frame`

- coupon_var:

  the variable for coupon identification

- x_var:

  the independent variable

- y_var:

  the dependent variable

- fn:

  a function defining the relationship between `Y` and `X`. See Details
  for more information.

- par:

  the initial guess for the parameters

- n_bins:

  the number of bins to average the data inside into before fitting

- method:

  The method to be used by
  [`optim()`](https://rdrr.io/r/stats/optim.html). Defaults to
  "L-BFGS-B"

- ...:

  extra parameters to be passed to
  [`optim()`](https://rdrr.io/r/stats/optim.html)

## Value

an object of class `average_curve_optim` with the following content:

- `data` the original data provided to the function

- `binned_data` the data after the binning/averaging operation

- `fn` the function supplied

- `fit_optim` the results of the call to `optim`

- `call` the call

- `n_bins` the number of bins specified by the user

- `max_x` the upper end of the range used for fitting

- `y_var` the independent (`Y`) variable

- `x_var` the dependent (`X`) variable

## Details

The function `fn` must have two arguments. The first argument must be
the value of the independent variable (`X`): this must be a numeric
value (of length one). The second argument must be a vector of the
parameters of the model, which are to be varied in order to obtain the
best fit. See below for an example.

## See also

[`optim()`](https://rdrr.io/r/stats/optim.html),
[`average_curve_lm()`](https://cmstatrExt.cmstatr.net/reference/average_curve_lm.md),
[`print.average_curve_optim()`](https://cmstatrExt.cmstatr.net/reference/print.average_curve_optim.md),
[`augment.average_curve_optim()`](https://cmstatrExt.cmstatr.net/reference/augment.average_curve_optim.md)

## Examples

``` r
# using the `pa12_tension` dataset and fitting a cubic polynomial with
# zero intercept:
curve_fit <- average_curve_optim(
  pa12_tension,
  Coupon,
  Strain,
  Stress,
  function(strain, par) {
    sum(par * c(strain, strain^2, strain^3))
  },
  c(c1 = 1, c2 = 1, c3 = 1),
  n_bins = 100
)
## Range: ` Strain ` in  [ 0,  0.1409409 ]
##
## Call:
## average_curve_optim(data = pa12_tension, coupon_var = Coupon,
##                     x_var = Strain, y_var = Stress,
##                     fn = function(strain, par) {
##                       sum(par * c(strain, strain^2, strain^3))
##                     }, par = c(c1 = 1, c2 = 1, c3 = 1), n_bins = 100)
##
## Parameters:
##       c1        c2        c3
## 1174.372 -8783.106 20585.898
```
