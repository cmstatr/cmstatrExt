#include <Rcpp.h>
// #include <cmath>
// #include "root.h"
// #include "integration.h"
#include "acceptance.h"
#include <testthat.h>
#include "testthat-exp.h"


//' p-Value for two-sample equivalency
//'
//' @description
//' Calculates the p-Value for a two-sample acceptance test.
//' This test considers the sample size of the qualification
//' sample (`n`) and the acceptance sample (`m`).
//'
//' Two test statistics are required:
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
//' @param n the size of the qualification sample
//' @param m the size of the acceptance sample
//' @param t1 the test statistic described above. May be a vector.
//' @param t2 the test statistic described above. May be a vector.
//'
//' @return a vector of p-Values of the same length as t1 and t2
//'
//' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector p_equiv_two_sample(int n, int m,
                                       Rcpp::NumericVector t1,
                                       Rcpp::NumericVector t2) {
  if (n < 3 || m < 3) {
    ::Rf_error("Both n and m must be 3 or greater");
  }
  
  if (t1.size() != t2.size()) {
    ::Rf_error("t1 and t2 must be of the same length");
  }
  const int num_vals = t1.size();
  for(int i = 0; i < num_vals; ++i) {
    if(t1[i] < t2[i]) {
      ::Rf_error("t2 must be less than t1");
    }
  }
  AcceptanceTwoSample an = AcceptanceTwoSample(n, m);
  Rcpp::NumericVector result = Rcpp::NumericVector(num_vals);
  
  for(int i = 0; i < num_vals; ++i) {
    result[i] = an.calc_p_value(t1[i], t2[i]);
  }
  
  return result;
}
