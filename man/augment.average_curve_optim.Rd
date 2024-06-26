% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/curve-fit.R
\name{augment.average_curve_optim}
\alias{augment.average_curve_optim}
\title{Augment a \code{data.frame} with the results from \code{average_curve_optim}}
\usage{
\method{augment}{average_curve_optim}(x, newdata = NULL, extrapolate = FALSE, ...)
}
\arguments{
\item{x}{an \code{average_curve_optim} object}

\item{newdata}{(optional) a new \code{data.frame} to which to augment the object}

\item{extrapolate}{whether to show the curve fit on all data or only
the data within the original fitted range. Default: FALSE}

\item{...}{ignored}
}
\value{
a \code{data.frame} with new columns \code{.fit}, \code{.extrapolate} and
\code{.residual}
}
\description{
Augment a \code{data.frame} with the results from \code{average_curve_optim}
}
\examples{
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
augment(curve_fit)
## # A tibble: 3,105 × 6
## Coupon     Strain  Stress  .fit .extrapolate .residual
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
## # i 3,095 more rows
## # i Use `print(n = ...)` to see more rows

}
\seealso{
\code{\link[=average_curve_lm]{average_curve_lm()}}
}
