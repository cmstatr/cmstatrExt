#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include <stdio.h>
#include <emscripten/emscripten.h>
#include "../root.h"
#include "../acceptance.h"
#include "power_sim.h"
#include "Rf_error.h"


Catch::Session session; // There must be exactly one instance

#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN
#endif
EXTERN EMSCRIPTEN_KEEPALIVE int do_tests() {
  int numFailed = session.run();
  return numFailed;
}

int main() {
  // Do nothing in main()
  // Various other function calls provide functionality
  return 0;
}


#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN
#endif
EXTERN EMSCRIPTEN_KEEPALIVE int k_equiv_two_sample(int n, int m, double alpha,
                                                   double* factors) {
  if (n < 3 || m < 3) {
    Rf_error("Both n and m must be 3 or greater (are %i and %i)", n, m);
    return -1;
  }
  
  AcceptanceTwoSample an = AcceptanceTwoSample(n, m);
  an.calculate_factors(alpha);
  factors[0] = an.k1;
  factors[1] = an.k2;
  return 0;
}


std::vector<double> range(double min, double max, size_t N) {
  std::vector<double> range;
  double delta = (max - min) / double(N - 1);
  for(int i=0; i<N; i++) {
    range.push_back(min + i * delta);
  }
  return range;
}


#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN
#endif
EXTERN EMSCRIPTEN_KEEPALIVE int power_mean(int n, int m, double k1, double k2,
                                           int steps, double* mu_pop_equiv, double* rejection_rate) {
  std::vector<double> sd_pop_equiv = range(1, 1, steps);
  std::vector<double> mu_pop_equiv_vec(steps);
  for(int i = 0; i < steps; ++i) {
    mu_pop_equiv_vec[i] = mu_pop_equiv[i];
  }
  std::vector<double> rej_rate_vec(steps);
  bool result = power_sim_dual_normal(
    n, m, 2500,
    0., 1.,
    mu_pop_equiv_vec, sd_pop_equiv,
    k1, k2, 
    rej_rate_vec
  );
  for(int i = 0; i < steps; ++i) {
    rejection_rate[i] = rej_rate_vec[i];
  }
  return result;
}

#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN
#endif
EXTERN EMSCRIPTEN_KEEPALIVE double p_equiv_two_sample(int n, int m,
                                                      double t1, double t2) {
  if (n < 3 || m < 3) {
    Rf_error("Both n and m must be 3 or greater (are %i and %i)", n, m);
    return -1.;
  }
  
  if(t1 < t2) {
    Rf_error("t2 must be less than t1 (are %f and %f)", t1, t2);
    return -1.;
  }
  AcceptanceTwoSample an = AcceptanceTwoSample(n, m);
  return an.calc_p_value(t1, t2);
}
