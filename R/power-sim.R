
#' Rejection rate for dual acceptance criteria based via simulation
#'
#' @description
#' Performs Monte Carlo simulation to determine the rejection rate of a
#' dual acceptance criteria (based on sample minimum and mean). By specifying
#' several sets of parameters for the "equivalency" distribution, a power
#' curve for the acceptance test can be determined.
#'
#' @param n_qual the sample size of the qualification sample
#' @param m_equiv the sample size of the equivalency/acceptance sample
#' @param replicates the number of simulated qualification samples and
#'                   equivalency samples. If a single value is given, the
#'                   same numbers used for both, if a vector of length two
#'                   is given, the first element is the number of qualification
#'                   replicates and the second element is the number of
#'                   equivalency replicates.
#' @param distribution a function name for generating a random sample
#'                     (defaults to "rnorm")
#' @param param_qual a data.frame (must be single row) with columns matching the
#'                   arguments of the distribution function
#' @param param_equiv a data.frame with columns matching the arguments of the
#'                    distribution function. The simulation is repeated with
#'                    the parameters specified in each row of the data.frame.
#' @param k1 a factor for determining the acceptance criterion for sample
#'           minimum, which is calculated as `mean_qual - k1 * sd_qual`
#' @param k2 a factor for determining the acceptance criterion for sample
#'           average, which is calculated as `mean_qual - k2 * sd_qual`
#'
#' @returns
#' A `data.frame`. The first column(s) are duplicate of the `data.frame`
#' passed in the argument `param_equiv`. The last column is named
#' `Rejection Rate` and has a value equal to the number of samples rejected
#' for each simulation run.
#'
#' @details
#' This function performs simulation to estimate the performance of the
#' dual acceptance criteria commonly used for composite materials in
#' aerospace applications. These criteria are based on setting lower limits
#' on the minimum individual (lower extremum) and the mean of an "acceptance"
#' sample. These limits are computed based on the sample mean and sample
#' standard deviation of an initial "qualification" sample. The criteria
#' are intended to be a test of non-inferiority to determine if the material
#' lots from which the "acceptance" samples are drawn should be accepted for
#' production. Another common use of these criteria are to determine if a
#' process change, equipment change, or second manufacturing site is
#' acceptable for production.
#'
#' For each set of distribution parameters given by the rows of `param_equiv`,
#' a number of samples of size `m_equiv` are generated using the function
#' `distribution`. Next, a number of qualification samples of size `n_qual`
#' are generated
#' using the `distribution` function and the parameters given in `param_qual`.
#' Limits for minimum and average are determined for each qualification
#' sample. Each equivalency sample is compared with the limits determined
#' from each qualification sample. The number of replicate in this simulation
#' is given by `replicates`: if this is a vector of length two, the first
#' element is the number of qualification samples and the second is the number
#' of equivalency samples; if `replicates` is a single value, the same number
#' is used for the number of qualification samples and acceptance samples.
#' Therefore, for each row of `param_equiv` a total of
#' `replicates[1] * replicates[2]` criteria are evaluated.
#'
#' The argument `distribution` must correspond with a function that generates
#' (pseudo) random numbers. This function must have an argument `n` that
#' specifies the sample size to be generated. When the argument `distribution`
#' matches certain common distribution functions (such as `rnorm`), the C++
#' implementation of the random number generation function is used instead of
#' calling R code, which results in a significant speedup.
#'
#' @examples
#' # Compute a power curve for a dual acceptance criteria for a qualification
#' # sample size of 18 and an equivalency sample size of 6, using 2000
#' # replicates. A standard normal distribution is used and the power to
#' # detect a decrease in mean is determined.
#' set.seed(12345) # make it reproducible
#' power_sim_dual(
#'   18, 6,
#'   2000,
#'   "rnorm",
#'   data.frame(mean = 0, sd = 1),
#'   data.frame(mean = c(-2, -1.5, -1, 0.5, 0), sd = 1),
#'   2.959, 0.954
#' )
#' ##   mean sd Rejection Rate
#' ## 1 -2.0  1     0.98349975
#' ## 2 -1.5  1     0.88186900
#' ## 3 -1.0  1     0.56382425
#' ## 4  0.5  1     0.00864025
#' ## 5  0.0  1     0.04826250
#'
#' @seealso [`k_equiv_two_sample`]
#'
#' @export
power_sim_dual <- function(n_qual, m_equiv,
                           replicates,
                           distribution = "rnorm",
                           param_qual,
                           param_equiv,
                           k1, k2) {
  stopifnot("param_qual must be a data.frame" = is.data.frame(param_qual))
  stopifnot("param_qual must have exactly one row" = nrow(param_qual) == 1)
  stopifnot("param_qual must have at least one column" = ncol(param_qual) > 0)
  stopifnot("param_equiv must be a data.frame" = is.data.frame(param_equiv))
  stopifnot("param_equiv must have at least one row" = nrow(param_equiv) >= 1)
  stopifnot("param_equiv must have at least one column" = ncol(param_equiv) > 0)

  known_distributions <- c(
    "norm", "rnorm"
  )

  if (as.character(distribution) %in% known_distributions) {
    return(power_sim_dual_generic(
      n_qual, m_equiv,
      replicates,
      distribution,
      function() NULL,
      param_qual,
      param_equiv,
      k1, k2
    ))
  } else {
    env <- parent.frame(n = 1)
    dist_fcn <- function(args, ii) {
      # The C++ code uses a zero-based index, R uses a one-based index, so we
      # need to correct that
      do.call(distribution, args[ii + 1, ], envir = env)
    }
    param_qual$n <- n_qual
    param_equiv$n <- m_equiv
    df <- power_sim_dual_generic(
      n_qual, m_equiv,
      replicates,
      "",
      dist_fcn,
      param_qual,
      param_equiv,
      k1, k2
    )
    # To be consistent with the built-in distributions, don't include the "n"
    # column
    df[["n"]] <- NULL
    return(df)
  }
}
