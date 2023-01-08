#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// TODO: allow for multiple methods
// TODO: allow testing various distributions
/*
 * Ideas:
 * - rep_qual and rep_equiv should be a vector of length 1 or 2
 * - mu_qual and sigma_qual should be a data.frame or named list
 * - mu_equiv and sigma_equiv should be a data.frame with one row per
 *   "delta"
 * - take parameter for distribution that can be matched to "norm",
 *   "weibull", etc., or an R function name.
 * 
 */

//' @export
// [[Rcpp::export]]
DataFrame power_sim(const int n_qual, const int m_equiv,
                    const int rep_qual, const int rep_equiv,
                    const double mu_qual, const double sigma_qual,
                    NumericVector mu_equiv, NumericVector sigma_equiv,
                    const double k1, const double k2) {
  if(mu_equiv.length() != sigma_equiv.length()) {
    ::Rf_error("mu_equiv and sigma_equiv must have the same length");
  }
  if(n_qual <= 0 || m_equiv <= 0) {
    ::Rf_error("n_qual and m_equiv must both be at least 1");
  }
  if(rep_qual <= 0 || rep_equiv <= 0) {
    ::Rf_error("rep_qual and rep_equiv must both be at least 1");
  }
  
  NumericVector min_equiv = NumericVector(rep_equiv);
  NumericVector avg_equiv = NumericVector(rep_equiv);
  IntegerVector accept_count = IntegerVector(mu_equiv.length());
  NumericVector reject_rate = NumericVector(mu_equiv.length());
  
  for(int j_param = 0; j_param < mu_equiv.length(); ++j_param) {
    const double mu = mu_equiv[j_param];
    const double sigma = sigma_equiv[j_param];
    
    for(int j_equiv = 0; j_equiv < rep_equiv; ++j_equiv) {
      const NumericVector x_equiv = ::rnorm(m_equiv, mu, sigma);
      min_equiv[j_equiv] = ::min(x_equiv);
      avg_equiv[j_equiv] = ::mean(x_equiv);
    }
    
    for(int j_qual = 0; j_qual < rep_qual; ++j_qual) {
      const NumericVector x_qual = ::rnorm(n_qual, mu_qual, sigma_qual);
      const double avg_qual = ::mean(x_qual);
      const double sd_qual = ::sd(x_qual);
      
      for(int j_equiv = 0; j_equiv < rep_equiv; ++j_equiv) {
        if(min_equiv[j_equiv] > avg_qual - k1 * sd_qual &&
           avg_equiv[j_equiv] > avg_qual - k2 * sd_qual) {
          accept_count[j_param]++;
        }
      }
    }
    
    reject_rate[j_param] = (double)(rep_qual * rep_equiv - accept_count[j_param]) /
      (rep_qual * rep_equiv);
  }
  
  return Rcpp::DataFrame::create(
    Named("mu_equiv") = mu_equiv,
    Named("sigma_equiv") = sigma_equiv,
    Named("Rejection Rate") = reject_rate
  );
}

