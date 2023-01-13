
#' @export
power_sim_dual <- function(n_qual, m_equiv,
                           replicates,
                           distribution,
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
