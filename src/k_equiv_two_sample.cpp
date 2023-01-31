#include <Rcpp.h>
// #include <cmath>
// #include "root.h"
// #include "integration.h"
#include "acceptance.h"
#include <testthat.h>
#include "testthat-exp.h"

//' Calculate the factors for a two-sample acceptance test
//'
//' @description
//' Calculates the factors k1 and k2, which are used for setting acceptance
//' values for lot acceptance. These factors consider both the
//' size of the qualification sample (`n`)
//' and the size of acceptance sample (`m`).
//' This test is detailed in a forthcoming paper.
//'
//' @param alpha the desired probability of Type 1 error
//' @param n the size of the qualification sample
//' @param m the size of the acceptance sample
//'
//' @return
//' A vector of length 2 with the contents `c(k1, k2)`
//' 
//' @references
//' Kloppenborg, S. (2023). Lot acceptance testing using sample mean and
//' extremum with finite qualification samples. Journal of Quality Technology,
//' https://doi.org/10.1080/00224065.2022.2147884
//'
//' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector k_equiv_two_sample(double alpha, int n, int m) {
  if (n < 3 || m < 3) {
    ::Rf_error("Both n and m must be 3 or greater");
  }
  
  AcceptanceTwoSample an = AcceptanceTwoSample(n, m);
  an.calculate_factors(alpha);
  
  return Rcpp::NumericVector::create(
    an.k1,
    an.k2
  );
}
