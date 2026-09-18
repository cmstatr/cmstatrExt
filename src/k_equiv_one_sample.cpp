#include <Rcpp.h>

#define _Rf_error Rcpp::stop
#define _Rf_warn Rcpp::warning

// #include <cmath>
// #include "root.h"
// #include "integration.h"
#include "acceptance.h"
#include <testthat.h>
#include "testthat-exp.h"

//' Calculate the factors for a one-sample acceptance test
//'
//' @description
//' Calculates the factors k1 and k2, which are used for setting acceptance
//' values for lot acceptance. These factors consider 
//' the size of acceptance sample (`m`).
//'
//' @param alpha the nominal significance of the test
//' @param m the size of the acceptance sample
//'
//' @return
//' A vector of length 2 with the contents `c(k1, k2)`
//' 
//' @details
//' This function is equivalent to [cmstatr::k_equiv()], but is implemented
//' in C++ instead of in R, and is hence slightly faster.
//' 
//' @references
//' Vangel, M. (2002). Lot Acceptance and Compliance Testing Using the
//' Sample Mean and an Extremum.
//' Technometrics. https://doi.org/10.1198/004017002188618428
//'
//' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector k_equiv_one_sample(double alpha, int m) {
  if (m < 3) {
    _Rf_error("m must be 3 or greater");
  }
  if (alpha <= 0.) {
    _Rf_error("alpha must be positive");
  }
  if (alpha >= 1.) {
    _Rf_error("alpha must be less than 1");
  }
  if (alpha < 1e-5 || alpha > 0.5) {
    _Rf_warn(
      "k-factor solution has only been validated for 1e-5 <= alpha <= 0.5");
  }
  
  AcceptanceVangel an = AcceptanceVangel(m);
  an.calculate_factors(alpha);
  
  return Rcpp::NumericVector::create(
    an.k1,
    an.k2
  );
}
