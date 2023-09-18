

#' Generate an average stress-strain (or bearing load-deformation) curve
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
#' @returns an object of class `average_curve` with the following content:
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
#' curve_fit <- average_curve(
#'   pa12_tension,
#'   Sample,
#'   Stress ~ I(Strain) + I(Strain^2) + I(Strain^3) + 0,
#'   n_bins = 100
#' )
#' print(curve_fit)
#' ## Range: ` Strain ` in  [ 0,  0.1409409 ]
#' ##
#' ## Call:
#' ##   average_curve(data = pa12_tension, coupon_var = Sample, model = Stress ~
#' ##                 I(Strain) + I(Strain^2) + I(Strain^3) + 0, n_bins = 100)
#' ##
#' ## Coefficients:
#' ##    I(Strain)   I(Strain^2)   I(Strain^3)
#' ##        1174         -8783         20586
#'
#' @seealso [`~`][base::~], [`I()`][base::I()]
#'
#' @importFrom dplyr mutate if_else summarise ungroup group_by select n filter
#' @importFrom dplyr n_groups
#' @importFrom stats lm
#' @importFrom rlang `:=` f_lhs f_rhs warn .data
#' @export
average_curve <- function(data, coupon_var, model, n_bins = 100) {
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

  fit <- lm(model, data = binned_dat)
  fit$call <- match.call()

  res <- list(
    data = data,
    binned_data = binned_dat,
    fit_lm = fit,
    n_bins = n_bins,
    max_x = max_x,
    y_var = y_var,
    x_var = x_var
  )
  class(res) <- "average_curve"
  res
}

#' Print an `average_curve` object
#'
#' @param x an `average_curve` object
#' @param ... additional arguments passed to `print.lm`
#'
#' @returns The object passed to this method is invisibly returned
#'          (via `invisible(x)`).
#'
#' @seealso [average_curve()]
#'
#' @method print average_curve
#' @export
print.average_curve <- function(x, ...) {
  cat("\nRange: `", x$x_var, "` in ",
      "[ 0, ", x$max_x, "]\n")
  print(x$fit_lm, ...)
  invisible(x)
}

#' Summarize an `average_curve` object
#'
#' @param object an `average_curve` object
#' @param ... arguments passed to `summary.lm`
#'
#' @returns No return value. This method only produces printed output.
#'
#' @seealso [average_curve()]
#'
#' @method summary average_curve
#' @export
summary.average_curve <- function(object, ...) {
  cat("\nRange: `", object$x_var, "` in ",
      "[0, ", object$max_x, "]\n")
  cat("n_bins = ", object$n_bins, "\n")
  summary(object$fit_lm, ...)
}

#' Augment a `data.frame` with the results from `average_curve`
#'
#' @param x an `average_curve` object
#' @param newdata (optional) a new `data.frame` to which to augment the object
#' @param extrapolate whether to show the curve fit on all data or only
#'                    the data within the original fitted range. Default: FALSE
#' @param ... ignored
#'
#' @returns a `data.frame` with new columns `.fit`, `.extrapolate` and
#'          `.residual`
#'
#' @examples
#' curve_fit <- average_curve(
#'   pa12_tension,
#'   Sample,
#'   Stress ~ I(Strain) + I(Strain^2) + I(Strain^3) + 0,
#'   n_bins = 100
#' )
#' augment(curve_fit)
#' ## # A tibble: 3,105 × 6
#' ##    Sample     Strain  Stress  .fit .extrapolate .residual
#' ##    <chr>       <dbl>   <dbl> <dbl> <lgl>            <dbl>
#' ##  1 Sample 4 0        -0.353  0     FALSE          -0.353
#' ##  2 Sample 4 0.000200 -0.0604 0.235 FALSE          -0.295
#' ##  3 Sample 4 0.000400  0.283  0.469 FALSE          -0.185
#' ##  4 Sample 4 0.000601  0.475  0.702 FALSE          -0.228
#' ##  5 Sample 4 0.000801  0.737  0.935 FALSE          -0.198
#' ##  6 Sample 4 0.00100   0.803  1.17  FALSE          -0.364
#' ##  7 Sample 4 0.00120   1.25   1.40  FALSE          -0.151
#' ##  8 Sample 4 0.00140   1.32   1.63  FALSE          -0.305
#' ##  9 Sample 4 0.00160   1.53   1.86  FALSE          -0.325
#' ## 10 Sample 4 0.00180   2.01   2.09  FALSE          -0.0735
#' ## # ℹ 3,095 more row
#' ## # ℹ Use `print(n = ...)` to see more rows
#'
#' @seealso [average_curve()]
#'
#' @method augment average_curve
#' @importFrom stats predict.lm
#' @export
augment.average_curve <- function(x, newdata = NULL, extrapolate = FALSE,
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
