#include <stdio.h>
#include <vector>
#include <algorithm>
#include <math.h>
#include "nmath/nmath.h"


std::vector<double> rnorm3(const int n, const double mu, const double sigma) {
  std::vector<double> result(n);
  for (int i = 0; i < n; i++)
  {
    result[i] = R::rnorm(mu, sigma);
  }
  return result;
}

void min_avg(std::vector<double> x, double &min, double &avg) {
  min = INFINITY;
  avg = 0.;
  for(double xi : x) {
    avg += xi;
    if(xi < min) {
      min = xi;
    }
  }
  avg /= x.size();
}

void avg_sd(std::vector<double> x, double &avg, double &sd) {
  avg = 0.;
  sd = 0.;
  for(double xi : x) {
    avg += xi;
  }
  avg /= x.size();

  for(double xi : x) {
    sd += (xi - avg) * (xi - avg);
  }
  sd /= x.size() - 1.;
  sd = sqrt(sd);
}


bool power_sim_dual_normal(
    const int n_qual, const int m_equiv,
    const int replicates,
    const double mu_pop_qual, const double sd_pop_qual,
    std::vector<double> mu_pop_equiv, std::vector<double> sd_pop_equiv,
    const double k1, const double k2,
    std::vector<double> &reject_rate) {
  
  if(n_qual <= 0 || m_equiv <= 0) {
    printf("n_qual and m_equiv must both be at least 1\n");
    return false;
  }
  int rep_qual = 0;
  int rep_equiv = 0;
  if(replicates <= 0) {
    printf("Number of replicates must be greater than zero\n");
    return false;
  }
  rep_qual = replicates;
  rep_equiv = replicates;
  if(mu_pop_equiv.size() != sd_pop_equiv.size()) {
    printf("mu_pop_equiv and sd_pop_equiv must be the same length\n");
    return false;
  }
  if(reject_rate.size() != mu_pop_equiv.size()) {
    printf("reject_rate is the wrong size.\n");
    return false;
  }
  
  std::vector<double> min_equiv(rep_equiv);
  std::vector<double> avg_equiv(rep_equiv);
  const int step_count = mu_pop_equiv.size();
  std::vector<int> accept_count(step_count);
  
  for(int j_param = 0; j_param < step_count; ++j_param) {
    accept_count[j_param] = 0;

    for(int j_equiv = 0; j_equiv < rep_equiv; ++j_equiv) {
      const std::vector<double> x_equiv = rnorm3(m_equiv, mu_pop_equiv[j_param], sd_pop_equiv[j_param]);
      min_avg(x_equiv, min_equiv[j_equiv], avg_equiv[j_equiv]);
    }
    
    for(int j_qual = 0; j_qual < rep_qual; ++j_qual) {
      const std::vector<double> x_qual = rnorm3(n_qual, mu_pop_qual, sd_pop_qual);
      double avg_qual;
      double sd_qual;
      avg_sd(x_qual, avg_qual, sd_qual);
      
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
  
  return true;
}
