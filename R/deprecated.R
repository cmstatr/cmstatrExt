#' p-Value for one-sample equivalency
#' 
#' @description
#' `r lifecycle::badge("deprecated")`
#' 
#' Use [p_equiv_one_sample()] instead of `p_equiv()`.
#' 
#' @param m the size of the acceptance sample
#' @param t1 the test statistic described above. May be a vector.
#' @param t2 the test statistic described above. May be a vector.
#' 
#' @importFrom lifecycle deprecate_warn
#' 
#' @export
p_equiv <- function(m, t1, t2) {
  deprecate_warn("0.5.0", "p_equiv()", "p_equiv_one_sample()")
  p_equiv_one_sample(m, t1, t2)
}
