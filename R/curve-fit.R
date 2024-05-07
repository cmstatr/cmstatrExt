

#' Generate an average curve using `lm`
#'
#' @description
#' The user must decide on a single dependent variable (`Y`) and a
#' single independent variable (`X`). The user will specify a `formula` with
#' the relationship between the dependent and independent variables.
#' For a `data.frame` containing stress-strain (or load-deflection) data for
#' more than one coupon, the maximum value of `X` for each coupon is found and
#' the smallest maximum value determines the range over which the curve
#' fit is performed: the range is from zero to this value. Only positive
#' values of `X` are considered. For each coupon individually, the data is
#' divided into a user-specified number of bins and averaged within each bin.
#' The resulting binned/averaged data is then passed to [stats::lm()] to perform
#' the curve fitting.
#'
#' @param data a `data.frame`
#' @param coupon_var the variable for coupon identification
#' @param model a `formula` for the curve to fit
#' @param n_bins the number of bins to average the data inside into before
#'               fitting
#'
#' @details
#' When specifying the formula (argument `model`), there are two things to
#' keep in mind. First, based on physical behavior, it is normally desirable
#' to set the intercept to zero (e.g. so that there is 0 stress at 0 strain).
#' To do this, include a term `+0` in the formula. Second, when specifying
#' a term for a power of the `X` variable (for example, $X^2$), this needs
#' to be wrapped inside the "as-is" operator `I()`, otherwise, `R` will
#' treat it as an interaction term, rather than an exponent. In other words,
#' if you want to include a quadratic term, you need to write `I(X^2)`
#' (replacing `X` with the appropriate variable from your `data.frame`).
#'
#' @returns an object of class `average_curve_lm` with the following content:
#' - `data` the original data provided to the function
#' - `binned_data` the data after the binning/averaging operation
#' - `fit_lm` the results of the call to `lm`
#' - `n_bins` the number of bins specified by the user
#' - `max_x` the upper end of the range used for fitting
#' - `y_var` the independent (`Y`) variable
#' - `x_var` the dependent (`X`) variable
#'
#' @examples
#' # using the `pa12_tension` dataset and fitting a cubic polynomial with
#' # zero intercept:
#' curve_fit <- average_curve_lm(
#'   pa12_tension,
#'   Coupon,
#'   Stress ~ I(Strain) + I(Strain^2) + I(Strain^3) + 0,
#'   n_bins = 100
#' )
#' print(curve_fit)
#' ## Range: ` Strain ` in  [ 0,  0.1409409 ]
#' ##
#' ## Call:
#' ##   average_curve_lm(data = pa12_tension, coupon_var = Coupon,
#' ##                    model = Stress ~ I(Strain) + I(Strain^2) + I(Strain^3)
#' ##                    + 0, n_bins = 100)
#' ##
#' ## Coefficients:
#' ##    I(Strain)   I(Strain^2)   I(Strain^3)
#' ##        1174         -8783         20586
#'
#' @seealso [`~`][base::~], [`I()`][base::I()], [`lm()`][stats::lm()],
#'          [average_curve_optim()], [print.average_curve_lm()],
#'          [summary.average_curve_lm()], [augment.average_curve_lm()]
#'
#' @importFrom dplyr mutate if_else summarise ungroup group_by select n filter
#' @importFrom dplyr n_groups
#' @importFrom stats lm
#' @importFrom rlang `:=` f_lhs f_rhs warn .data
#' @export
average_curve_lm <- function(data, coupon_var, model, n_bins = 100) {
  if (length(all.vars(f_lhs(model))) != 1) {
    stop("LHS of `model` must have exactly one variable")
  }
  if (length(all.vars(f_rhs(model))) != 1) {
    stop("RHS of `model` must have exactly one variable")
  }
  if (all.vars(f_lhs(model)) == all.vars(f_rhs(model))) {
    stop("RHS and LHS of `model` must use different variables")
  }

  y_var <- as.name(all.vars(f_lhs(model)))
  x_var <- as.name(all.vars(f_rhs(model)))

  binned_dat <- average_curve_bin_data(data, {{ coupon_var }},
                                       {{ x_var }}, {{ y_var }}, n_bins)

  fit <- lm(model, data = binned_dat$binned_dat)
  fit$call <- match.call()

  res <- list(
    data = data,
    binned_data = binned_dat$binned_dat,
    fit_lm = fit,
    n_bins = n_bins,
    max_x = binned_dat$max_x,
    y_var = y_var,
    x_var = x_var
  )
  class(res) <- "average_curve_lm"
  res
}

#' Generate an average curve using `optim`
#'
#' @description
#' The user must decide on a single dependent variable (`Y`) and a
#' single independent variable (`X`). The user will specify a function defining
#' the relationship between the dependent and independent variables.
#' For a `data.frame` containing stress-strain (or load-deflection) data for
#' more than one coupon, the maximum value of `X` for each coupon is found and
#' the smallest maximum value determines the range over which the curve
#' fit is performed: the range is from zero to this value. Only positive
#' values of `X` are considered. For each coupon individually, the data is
#' divided into a user-specified number of bins and averaged within each bin.
#' The resulting binned/averaged data is then used for curve fitting.
#' The mean squared error between the observed value of `Y` and the result of
#' the user-specified function evaluated at each `X` is minimized by varying
#' the parameters `par`.
#'
#' @param data a `data.frame`
#' @param coupon_var the variable for coupon identification
#' @param x_var the independent variable
#' @param y_var the dependent variable
#' @param fn a function defining the relationship between `Y` and `X`.
#'           See Details for more information.
#' @param par the initial guess for the parameters
#' @param n_bins the number of bins to average the data inside into before
#'               fitting
#' @param method The method to be used by `optim()`. Defaults to "L-BFGS-B"
#' @param ... extra parameters to be passed to `optim()`
#'
#' @details
#' The function `fn` must have two arguments. The first argument must be the
#' value of the independent variable (`X`): this must be a numeric value
#' (of length one). The second argument must be a vector of the parameters of
#' the model, which are to be varied in order to obtain the best fit. See below
#' for an example.
#'
#' @returns an object of class `average_curve_optim` with the following content:
#' - `data` the original data provided to the function
#' - `binned_data` the data after the binning/averaging operation
#' - `fn` the function supplied
#' - `fit_optim` the results of the call to `optim`
#' - `call` the call
#' - `n_bins` the number of bins specified by the user
#' - `max_x` the upper end of the range used for fitting
#' - `y_var` the independent (`Y`) variable
#' - `x_var` the dependent (`X`) variable
#'
#' @examples
#' # using the `pa12_tension` dataset and fitting a cubic polynomial with
#' # zero intercept:
#' curve_fit <- average_curve_optim(
#'   pa12_tension,
#'   Coupon,
#'   Strain,
#'   Stress,
#'   function(strain, par) {
#'     sum(par * c(strain, strain^2, strain^3))
#'   },
#'   c(c1 = 1, c2 = 1, c3 = 1),
#'   n_bins = 100
#' )
#' ## Range: ` Strain ` in  [ 0,  0.1409409 ]
#' ##
#' ## Call:
#' ## average_curve_optim(data = pa12_tension, coupon_var = Coupon,
#' ##                     x_var = Strain, y_var = Stress,
#' ##                     fn = function(strain, par) {
#' ##                       sum(par * c(strain, strain^2, strain^3))
#' ##                     }, par = c(c1 = 1, c2 = 1, c3 = 1), n_bins = 100)
#' ##
#' ## Parameters:
#' ##       c1        c2        c3
#' ## 1174.372 -8783.106 20585.898
#'
#' @seealso [`optim()`][stats::optim()], [average_curve_lm()],
#'          [print.average_curve_optim()], [augment.average_curve_optim()]
#'
#' @importFrom dplyr mutate if_else summarise ungroup group_by select n filter
#' @importFrom dplyr n_groups
#' @importFrom rlang `:=` f_lhs f_rhs warn .data enquo0 as_label
#' @importFrom stats optim
#' @export
average_curve_optim <- function(data, coupon_var, x_var, y_var,
                                fn, par, n_bins = 100, method = "L-BFGS-B",
                                ...) {

  test_fn_result <- fn(0, par)
  if (length(test_fn_result) != 1 || !is.numeric(test_fn_result)) {
    stop("The return type of `fn` must be `numeric` of length 1")
  }

  binned_dat <- average_curve_bin_data(data, {{ coupon_var }},
                                       {{ x_var }}, {{ y_var }}, n_bins)

  x_vec <- binned_dat$binned_dat[[as_label(enquo0(x_var))]]
  y_vec <- binned_dat$binned_dat[[as_label(enquo0(y_var))]]
  fn_ss <- function(par_i) {
    y_prime <- vapply(x_vec, function(x) fn(x, par_i), numeric(1))
    sum((y_prime - y_vec) ^ 2)
  }

  optim_res <- optim(par, fn_ss, method = method, ...)

  if (optim_res$convergence != 0) {
    warn(paste0("`optim` failed to converge: ", optim_res$message))
  }

  res <- list(
    data = data,
    binned_data = binned_dat$binned_dat,
    fn = fn,
    fit_optim = optim_res,
    call = match.call(),
    n_bins = n_bins,
    max_x = binned_dat$max_x,
    y_var = as.symbol(as_label(enquo0(y_var))),
    x_var = as.symbol(as_label(enquo0(x_var)))
  )
  class(res) <- "average_curve_optim"
  res
}

#' @importFrom dplyr mutate if_else summarise ungroup group_by select n filter
#' @importFrom dplyr n_groups
#' @importFrom rlang `:=` f_lhs f_rhs warn .data
average_curve_bin_data <- function(data, coupon_var, x_var, y_var, n_bins) {
  max_x <- group_by(data, {{coupon_var}})
  max_x <- summarise(max_x, "maxx" = max({{x_var}}), "minx" = min({{x_var}}))
  max_x <- select(max_x, c({{coupon_var}}, "maxx", "minx"))
  coupon_ids <- max_x[[1]]
  if (any(max_x$maxx <= 0)) {
    stop("No positive values of the `X` variable for one or more coupon.
         Negative values are ignored.")
  }
  if (any(max_x$minx < 0)) {
    warn("Negative values of the `X` variable are ignored.")
  }
  max_x <- min(max_x$maxx)

  binned_dat <- group_by(data, {{coupon_var}})
  binned_dat <- dplyr::filter(binned_dat, {{x_var}} <= max_x)
  binned_dat <- ungroup(binned_dat)
  binned_dat <- group_by(binned_dat, {{coupon_var}})
  binned_dat <- mutate(binned_dat,
                       `.bin` = cut({{x_var}}, seq(0,
                                                   max_x,
                                                   length.out = n_bins)))
  binned_dat <- group_by(binned_dat, .data[[".bin"]], {{coupon_var}})
  binned_dat <- summarise(binned_dat, {{y_var}} := mean({{y_var}}),
                          {{x_var}} := mean({{x_var}}),
                          ".n" = n(),
                          .groups = "drop")

  if (any(binned_dat$`.n` == 0)) {
    warn(paste0(
      "Some bins for coupon(s) ",
      unique(binned_dat[[coupon_var]][binned_dat$`.n` == 0]),
      " are empty. This may indicate too many bins, ",
      "or missing portions of data."
    ))
  } else if (nrow(binned_dat) != n_bins * length(coupon_ids)) {
    warn("Some bins are empty. This may indicate too many bins or missing data")
  }
  binned_dat <- select(binned_dat, -c(".bin", ".n"))

  list(
    binned_dat = binned_dat,
    max_x = max_x
  )
}

#' Print an `average_curve_lm` object
#'
#' @param x an `average_curve_lm` object
#' @param ... additional arguments passed to `print.lm`
#'
#' @returns The object passed to this method is invisibly returned
#'          (via `invisible(x)`).
#'
#' @seealso [average_curve_lm()]
#'
#' @method print average_curve_lm
#' @export
print.average_curve_lm <- function(x, ...) {
  cat("\nRange: `", x$x_var, "` in ",
      "[ 0, ", x$max_x, "]\n")
  print(x$fit_lm, ...)
  invisible(x)
}

#' Print an `average_curve_optim` object
#'
#' @param x an `average_curve_optim` object
#' @param ... not used
#'
#' @returns The object passed to this method is invisibly returned
#'          (via `invisible(x)`).
#'
#' @seealso [average_curve_optim()]
#'
#' @method print average_curve_optim
#' @export
print.average_curve_optim <- function(x, ...) {
  cat("\nRange: `", x$x_var, "` in ",
      "[ 0, ", x$max_x, "]\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  if (x$fit_optim$convergence != 0) {
    cat("Failed to converge: ", x$fit_optim$message, "\n")
  }
  cat("Parameters:\n")
  print(x$fit_optim$par)
  invisible(x)
}

#' Summarize an `average_curve_lm` object
#'
#' @param object an `average_curve_lm` object
#' @param ... arguments passed to `summary.lm`
#'
#' @returns No return value. This method only produces printed output.
#'
#' @seealso [average_curve_lm()]
#'
#' @method summary average_curve_lm
#' @export
summary.average_curve_lm <- function(object, ...) {
  cat("\nRange: `", object$x_var, "` in ",
      "[0, ", object$max_x, "]\n")
  cat("n_bins = ", object$n_bins, "\n")
  summary(object$fit_lm, ...)
}

#' Augment a `data.frame` with the results from `average_curve_lm`
#'
#' @param x an `average_curve_lm` object
#' @param newdata (optional) a new `data.frame` to which to augment the object
#' @param extrapolate whether to show the curve fit on all data or only
#'                    the data within the original fitted range. Default: FALSE
#' @param ... ignored
#'
#' @returns a `data.frame` with new columns `.fit`, `.extrapolate` and
#'          `.residual`
#'
#' @examples
#' curve_fit <- average_curve_lm(
#'   pa12_tension,
#'   Coupon,
#'   Stress ~ I(Strain) + I(Strain^2) + I(Strain^3) + 0,
#'   n_bins = 100
#' )
#' augment(curve_fit)
#' ## # A tibble: 3,105 × 6
#' ##    Coupon     Strain  Stress  .fit .extrapolate .residual
#' ##    <chr>       <dbl>   <dbl> <dbl> <lgl>            <dbl>
#' ##  1 Coupon 4 0        -0.353  0     FALSE          -0.353
#' ##  2 Coupon 4 0.000200 -0.0604 0.235 FALSE          -0.295
#' ##  3 Coupon 4 0.000400  0.283  0.469 FALSE          -0.185
#' ##  4 Coupon 4 0.000601  0.475  0.702 FALSE          -0.228
#' ##  5 Coupon 4 0.000801  0.737  0.935 FALSE          -0.198
#' ##  6 Coupon 4 0.00100   0.803  1.17  FALSE          -0.364
#' ##  7 Coupon 4 0.00120   1.25   1.40  FALSE          -0.151
#' ##  8 Coupon 4 0.00140   1.32   1.63  FALSE          -0.305
#' ##  9 Coupon 4 0.00160   1.53   1.86  FALSE          -0.325
#' ## 10 Coupon 4 0.00180   2.01   2.09  FALSE          -0.0735
#' ## # i 3,095 more row
#' ## # i Use `print(n = ...)` to see more rows
#'
#' @seealso [average_curve_lm()]
#'
#' @method augment average_curve_lm
#' @importFrom stats predict.lm
#' @export
augment.average_curve_lm <- function(x, newdata = NULL, extrapolate = FALSE,
                                     ...) {
  if (is.null(newdata)) {
    newdata <- ungroup(x$data)
  }

  result <- mutate(
    ungroup(newdata),
    ".fit" = predict.lm(x$fit_lm, newdata)
  )
  result <- mutate(
    result,
    `.extrapolate` = !!x$x_var > x$max_x | !!x$x_var < 0
  )
  if (!extrapolate) {
    result <- mutate(
      result,
      `.fit` = if_else(.data[[".extrapolate"]], NA, .data[[".fit"]])
    )
  }
  result <- mutate(
    result,
    `.residual` = !!x$y_var - .data[[".fit"]]
  )

  result
}

#' Augment a `data.frame` with the results from `average_curve_optim`
#'
#' @param x an `average_curve_optim` object
#' @param newdata (optional) a new `data.frame` to which to augment the object
#' @param extrapolate whether to show the curve fit on all data or only
#'                    the data within the original fitted range. Default: FALSE
#' @param ... ignored
#'
#' @returns a `data.frame` with new columns `.fit`, `.extrapolate` and
#'          `.residual`
#'
#' @examples
#' curve_fit <- average_curve_optim(
#'   pa12_tension,
#'   Coupon,
#'   Strain,
#'   Stress,
#'   function(strain, par) {
#'     sum(par * c(strain, strain^2, strain^3))
#'   },
#'   c(c1 = 1, c2 = 1, c3 = 1),
#'   n_bins = 100
#' )
#' augment(curve_fit)
#' ## # A tibble: 3,105 × 6
#' ## Coupon     Strain  Stress  .fit .extrapolate .residual
#' ##    <chr>       <dbl>   <dbl> <dbl> <lgl>            <dbl>
#' ##  1 Coupon 4 0        -0.353  0     FALSE          -0.353
#' ##  2 Coupon 4 0.000200 -0.0604 0.235 FALSE          -0.295
#' ##  3 Coupon 4 0.000400  0.283  0.469 FALSE          -0.185
#' ##  4 Coupon 4 0.000601  0.475  0.702 FALSE          -0.228
#' ##  5 Coupon 4 0.000801  0.737  0.935 FALSE          -0.198
#' ##  6 Coupon 4 0.00100   0.803  1.17  FALSE          -0.364
#' ##  7 Coupon 4 0.00120   1.25   1.40  FALSE          -0.151
#' ##  8 Coupon 4 0.00140   1.32   1.63  FALSE          -0.305
#' ##  9 Coupon 4 0.00160   1.53   1.86  FALSE          -0.325
#' ## 10 Coupon 4 0.00180   2.01   2.09  FALSE          -0.0735
#' ## # i 3,095 more rows
#' ## # i Use `print(n = ...)` to see more rows
#'
#' @seealso [average_curve_lm()]
#'
#' @method augment average_curve_optim
#' @export
augment.average_curve_optim <- function(x, newdata = NULL, extrapolate = FALSE,
                                        ...) {
  if (is.null(newdata)) {
    newdata <- ungroup(x$data)
  }
  result <- mutate(
    ungroup(newdata),
    ".fit" = vapply(!!x$x_var, function(xi) {
      x$fn(xi, x$fit_optim$par)
    },
    numeric(1)
    )
  )
  result <- mutate(
    result,
    `.extrapolate` = !!x$x_var > x$max_x | !!x$x_var < 0
  )
  if (!extrapolate) {
    result <- mutate(
      result,
      `.fit` = if_else(.data[[".extrapolate"]], NA, .data[[".fit"]])
    )
  }
  result <- mutate(
    result,
    `.residual` = !!x$y_var - .data[[".fit"]]
  )
  result
}
