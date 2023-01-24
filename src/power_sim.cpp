#include <Rcpp.h>
using namespace Rcpp;



void check_columns(DataFrame param, std::string name,
                   std::vector<std::string> columns) {
  for (const std::string& col : columns) {
    if(!param.containsElementNamed(col.c_str())) {
      ::Rf_error("%s must contain the column %s", name.c_str(), col.c_str());
    }
  }
}

// [[Rcpp::export]]
DataFrame power_sim_dual_generic(
    const int n_qual, const int m_equiv,
    NumericVector replicates,
    std::string distribution,
    Function dist_function,
    DataFrame param_qual,
    DataFrame param_equiv,
    const double k1, const double k2) {
  
  if(n_qual <= 0 || m_equiv <= 0) {
    ::Rf_error("n_qual and m_equiv must both be at least 1");
  }
  int rep_qual = 0;
  int rep_equiv = 0;
  if(replicates.length() == 1) {
    if(replicates[0] <= 0) {
      ::Rf_error("Number of replicates must be greater than zero");
    }
    rep_qual = replicates[0];
    rep_equiv = replicates[0];
  } else if(replicates.length() == 2) {
    if(replicates[0] <= 0 || replicates[1] <= 0) {
      ::Rf_error("Number of replicates must be greater than zero");
    }
    rep_qual = replicates[0];
    rep_equiv = replicates[1];
  } else {
    ::Rf_error("replicates must be a single number or a vector of length 2");
  }
  if(param_qual.rows() != 1) {
    ::Rf_error("`param_qual` must have exactly one row");
  }
  
  std::function<NumericVector(const int, const DataFrame, const int)>
    rgenerate = NULL;
  
  if(distribution == "norm" || distribution == "rnorm") {
    check_columns(param_qual, "param_qual", {"mean", "sd"});
    check_columns(param_equiv, "param_equiv", {"mean", "sd"});
    rgenerate = [](const int n, const DataFrame param, const int i) {
      const NumericVector m = param["mean"], s = param["sd"];
      return ::rnorm(n, m[i], s[i]);
    };
  } else if (distribution == "") {
    rgenerate = [dist_function](const int n, const DataFrame param, const int i) {
      return dist_function(param, i);
    };
  } else {  // TODO: add other distributions
    ::Rf_error("Unsupported distribution function");
  }
  
  NumericVector min_equiv = NumericVector(rep_equiv);
  NumericVector avg_equiv = NumericVector(rep_equiv);
  const int step_count = param_equiv.rows();
  IntegerVector accept_count = IntegerVector(step_count);
  NumericVector reject_rate = NumericVector(step_count);
  
  for(int j_param = 0; j_param < step_count; ++j_param) {
    
    for(int j_equiv = 0; j_equiv < rep_equiv; ++j_equiv) {
      const NumericVector x_equiv = rgenerate(m_equiv, param_equiv, j_param);
      min_equiv[j_equiv] = ::min(x_equiv);
      avg_equiv[j_equiv] = ::mean(x_equiv);
    }
    
    for(int j_qual = 0; j_qual < rep_qual; ++j_qual) {
      const NumericVector x_qual = rgenerate(n_qual, param_qual, 0);
      const double avg_qual = ::mean(x_qual);
      const double sd_qual = ::sd(x_qual);
      
      for(int j_equiv = 0; j_equiv < rep_equiv; ++j_equiv) {
        if(min_equiv[j_equiv] > avg_qual - k1 * sd_qual &&
           avg_equiv[j_equiv] > avg_qual - k2 * sd_qual) {
          accept_count[j_param]++;
        }
      }
    }
    
    reject_rate[j_param] = (double)
      (rep_qual * rep_equiv - accept_count[j_param]) /
      (rep_qual * rep_equiv);
  }
  
  DataFrame retVal = Rcpp::clone(param_equiv);
  retVal["Rejection Rate"] = reject_rate;
  retVal.attr("class") = "data.frame";
  retVal.attr("row.names") = seq(1, step_count);
  return retVal;
}

