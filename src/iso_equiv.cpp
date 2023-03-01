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

//' @export
// [[Rcpp::export]]
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

