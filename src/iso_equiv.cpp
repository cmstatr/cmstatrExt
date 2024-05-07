#include <Rcpp.h>
#include <vector>
#include <cmath>
#include "root.h"
#include "acceptance.h"
using namespace Rcpp;

#ifndef WASM
#include <Rcpp.h>
#include <testthat.h>

#else // WASM

#include "wasm/nmath/nmath.h"
#include "wasm/catch.hpp"
#include "wasm/testthat_catch.h"
#include "wasm/Rf_error.h"

#endif // WASM

#include "testthat-exp.h"


double t_interpolate(const double x, const double a, const double b) {
  return a + x * (b - a);
}

std::vector<double> open_range(double min, double max, size_t N,
                               bool include_start) {
  std::vector<double> range;
  double delta = (max - min) / double(N);
  if (!include_start) min += 0.5 * delta;
  for(size_t i=0; i<N; i++) {
    range.push_back(min + i * delta);
  }
  return range;
}

//' Calculate t1 and t2 pairs that have the same p-Value
//'
//' @description
//' Calculates pairs of t1 and t2 values, which have the same p-value for the
//' two-sample equivalency test. See [p_equiv_two_sample()].
//'
//' @param n the size of the qualification sample
//' @param m the size of the equivalency sample
//' @param alpha the desired p-value
//' @param t1max the maximum value of t1 (only approximate)
//' @param t2max the maximum value of t2 (only approximate)
//' @param n_points the number of returned points is twice n_points
//'
//' @details
//' The values t1 and t2 are based on the transformation:
//' 
//' t1 = (X_mean - Y_min) / S
//'
//' t2 = (X_mean - Y_mean) / S
//'
//' Where:
//' - X_mean is the mean of the qualification sample
//' - S is the standard deviation of the qualification sample
//' - Y_min is the minimum from the acceptance sample
//' - Y_mean is the mean of the acceptance sample
//'
//' @return
//' A `data.frame` with values of t1 and t2
//' 
//' @examples
//'\donttest{
//' if(requireNamespace("tidyverse")){
//'   library(cmstatrExt)
//'   library(tidyverse)
//'   curve <- iso_equiv_two_sample(24, 8, 0.05, 4, 1.5, 10)
//'   curve
//' 
//'   curve %>%
//'     ggplot(aes(x = t1, y = t2)) +
//'       geom_path() +
//'       ggtitle("Acceptance criteria for alpha=0.05")
//'}
//'}
//'
//' @seealso
//' [p_equiv_two_sample()], [k_equiv_two_sample()]
//' 
//' @references
//' Kloppenborg, S. (2023). Lot acceptance testing using sample mean and
//' extremum with finite qualification samples. Journal of Quality Technology,
//' https://doi.org/10.1080/00224065.2022.2147884
//'
//' @export
// [[Rcpp::export(rng = false)]]
DataFrame iso_equiv_two_sample(const int n, const int m, const double alpha,
                               double t1max, double t2max,
                               const double n_points) {
  AcceptanceTwoSample an = AcceptanceTwoSample(n, m);
  std::vector<double> result_t1;
  std::vector<double> result_t2;

  // t1 and t2 will be bounded by the asymptotic values.
  const double t1asy = an.calc_r1(alpha);
  const double t2asy = an.calc_r2(alpha);

  // make sure that t1max is at least a bit bigger than t1asy
  t1max = fmax(t1max, 1.1 * t1asy);
  // enforce t1 >= t2
  t2max = fmin(t2max, t1asy);

  for(double t1b : open_range(t1max, t1asy, n_points, true)) {
    const double t2b = t2asy;
    const auto f = [t1b, t2b, t1max, t2max, an, alpha](const double x) {
      return alpha - an.calc_p_value(
        t_interpolate(x, t1max, t1b),
        t_interpolate(x, t2max, t2b));
    };
    double x;
    int retval = bisection(f, 0., 1., &x);
    if(retval != ROOT_RESULT_SUCCESS) {
      ::Rf_error("Failed to find root");
    }
    result_t1.push_back(t_interpolate(x, t1max, t1b));
    result_t2.push_back(t_interpolate(x, t2max, t2b));
  }
  
  for(double t2b : open_range(t2asy, t2max, n_points, false)) {
    const double t1b = t1asy;
    const auto f = [t1b, t2b, t1max, t2max, an, alpha](const double x) {
      return alpha - an.calc_p_value(
          t_interpolate(x, t1max, t1b),
          t_interpolate(x, t2max, t2b));
    };
    double x;
    int retval = bisection(f, 0., 1., &x);
    if(retval != ROOT_RESULT_SUCCESS) {
      ::Rf_error("Failed to find root");
    }
    result_t1.push_back(t_interpolate(x, t1max, t1b));
    result_t2.push_back(t_interpolate(x, t2max, t2b));
  }
  
  DataFrame curve;
  curve["t1"] = result_t1;
  curve["t2"] = result_t2;
  curve.attr("class") = "data.frame";
  curve.attr("row.names") = seq(1, n_points * 2);
  return curve;
}

