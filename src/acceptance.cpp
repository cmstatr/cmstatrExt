#include <Rcpp.h>
#include <cmath>
#include "root.h"
#include "integration.h"
#include "acceptance.h"
#include "doctest.h"
#include "doctest-ext.h"


AcceptanceBase::AcceptanceBase(const double m) :
  m{m},
  a_int{} {
  a_int.init([this](const double t) { return a_fcn(t); }, true);
}

double AcceptanceBase::h(const double t) {
  return exp(DNORM(t, true) - PNORM(t, false, true));
}

double AcceptanceBase::a_fcn(const double t) {
  const double hmt = t > 60. ? pow(t, -1.) : h(t) - t;
  return exp(
    - (m - 1.) * (DNORM(t, true) - PNORM(t, false, true)) +
      pow(m - 1., 2.) / (2. * m) * pow(hmt, 2.) +
      (m - 1.) * t * hmt
  ) * sqrt(t > 60. ? pow(t, -2.) : 1. - h(t) * hmt);
}


double AcceptanceBase::calc_lambda(const double t1,
                                   const double t2, const double x0) {
  const auto f = [this, t1, t2](const double lam) {
    return (m - 1.) / m * (h(lam) - lam) - t2 + t1;
  };
  const auto dfdlam = [this](const double lam) {
    const double pn = PNORM(lam, false, false);
    return (m - 1.) / m * (
        pow(DNORM(lam, false) / pn, 2.)
      - lam * DNORM(lam, false) / pn
      - 1.
    );
  };
  
  double result = 0.;
  int retval;

  retval = root(f, dfdlam, x0, &result);
  
  if (retval != ROOT_RESULT_SUCCESS)
  {
    const int retval_bisection = bisection(f, -1000, 1000, &result, 1000);
    if (retval_bisection != ROOT_RESULT_SUCCESS) {
      ::Rf_error("Root failed. (Newton code=%i, bisection code=%i)",
                 retval, retval_bisection);
    }
  }
  
  return result;
}

double AcceptanceBase::calc_f_joint_vangel(const double t1, const double t2) {
  const auto a_m = [this](const double t) { return a_fcn(t); };
  const auto f1 = [this, t2](const double t) {
    return 1.;
  };
  const auto f2 = [this, t1](const double t) {
    return PNORM(
      sqrt(m) * (t1 + (m - 1.) / m * (h(t) - t)),
      true, false
    );
  };
  
  const double lam = calc_lambda(t1, t2, 0.);
  
  const IntegrationMultOneInf num1 = IntegrationMultOneInf(
    a_m, f1, &a_int, -1., lam);
  const IntegrationMultOneInf num2 = IntegrationMultOneInf(
    a_m, f2, &a_int, +1., lam);
  
  return (PNORM(sqrt(m) * t2, true, false) * num1.result + num2.result) /
    a_int.result;
}

double AcceptanceVangel::calc_f_min(const double t1) {
  return 1. - pow(PNORM(t1, false, false), m);
}

double AcceptanceVangel::calc_f_mean(const double t2) {
  return PNORM(sqrt(m) * t2, true, false);
}

AcceptanceVangel::AcceptanceVangel(const double m)
  : AcceptanceBase(m) {
}

void AcceptanceVangel::calculate_factors(const double alpha) {
  auto calc_t2 = [this](const double t1) {
    return -QNORM(
        1. - pow(PNORM(-t1, false, false), this->m),
        true, false
    ) / sqrt(this->m);
  };
  
  auto f = [this, calc_t2, alpha](const double t1) {
    double t2 = calc_t2(t1);
    const double fx1 = calc_f_min(-t1);
    const double fxbar = calc_f_mean(-t2);
    const double fjoint = calc_f_joint_vangel(-t1, -t2);
    return fx1 + fxbar - fjoint - alpha;
  };
  
  bisection(f, -0.1, 1, &k1, 100);
  k2 = calc_t2(k1);
}

double AcceptanceVangel::calc_p_value(const double r1, const double r2) {
  const double fx1 = calc_f_min(-r1);
  const double fxbar = calc_f_mean(-r2);
  const double fjoint = calc_f_joint_vangel(-r1, -r2);
  return fx1 + fxbar - fjoint;
}

TEST_CASE("Acceptance Vangel") {
  SUBCASE("m=5, alpha=0.05") {
    AcceptanceVangel ag = AcceptanceVangel(5);
    ag.calculate_factors(0.05);
    CHECK_ALMOST_EQ(ag.k1, 2.5286, 0.001);
    CHECK_ALMOST_EQ(ag.k2, 0.8525, 0.001);
    CHECK_ALMOST_EQ(ag.calc_p_value(2.5286, 0.8525), 0.05, 1e-5);
  }
  SUBCASE("m=10, alpha=0.05") {
    AcceptanceVangel ag = AcceptanceVangel(10);
    ag.calculate_factors(0.05);
    CHECK_ALMOST_EQ(ag.k1, 2.7772, 0.001);
    CHECK_ALMOST_EQ(ag.k2, 0.6089, 0.001);
    CHECK_ALMOST_EQ(ag.calc_p_value(2.7772, 0.6089), 0.05, 1e-5);
  }
  SUBCASE("m=5, alpha=0.5") {
    AcceptanceVangel ag = AcceptanceVangel(5);
    ag.calculate_factors(0.5);
    CHECK_ALMOST_EQ(ag.k1, 1.3498, 0.001);
    CHECK_ALMOST_EQ(ag.k2, 0.1473, 0.001);
    CHECK_ALMOST_EQ(ag.calc_p_value(1.3498, 0.1473), 0.5, 1e-5);
  }
  SUBCASE("m=10, alpha=0.5") {
    AcceptanceVangel ag = AcceptanceVangel(10);
    ag.calculate_factors(0.5);
    CHECK_ALMOST_EQ(ag.k1, 1.7258, 0.001);
    CHECK_ALMOST_EQ(ag.k2, 0.1217, 0.001);
    CHECK_ALMOST_EQ(ag.calc_p_value(1.7258, 0.1217), 0.5, 1e-4);
  }
  SUBCASE("m=5, alpha=0.0005") {
    AcceptanceVangel ag = AcceptanceVangel(5);
    ag.calculate_factors(0.0005);
    CHECK_ALMOST_EQ(ag.k1, 3.8864, 0.05);
    CHECK_ALMOST_EQ(ag.k2, 1.5546, 0.05);
    CHECK_ALMOST_EQ(ag.calc_p_value(3.8864, 1.5546), 0.0005, 1e-6);
  }
  SUBCASE("m=10, alpha=0.0005") {
    AcceptanceVangel ag = AcceptanceVangel(10);
    ag.calculate_factors(0.0005);
    CHECK_ALMOST_EQ(ag.k1, 4.0541, 0.05);
    CHECK_ALMOST_EQ(ag.k2, 1.1002, 0.05);
    CHECK_ALMOST_EQ(ag.calc_p_value(4.0541, 1.1002), 0.0005, 1e-6);
  }
}


//' p-Value for one-sample equivalency
//'
//' @description
//' Calculates the p-Value for a one-sample acceptance test
//' based on Vangel (2002).
//' This test considers the sample size of the acceptance sample (`m`).
//'
//' Two test statistics are required:
//'
//' t1 = (mu - Y_min) / sigma
//'
//' t2 = (mu - Y_mean) / sigma
//'
//' Where:
//' - mu is the mean of the population
//' - sigma is the standard deviation of the population
//' - Y_min is the minimum from the acceptance sample
//' - Y_mean is the mean of the acceptance sample
//'
//' @param m the size of the acceptance sample
//' @param t1 the test statistic described above. May be a vector.
//' @param t2 the test statistic described above. May be a vector.
//'
//' @return a vector of p-Values of the same length as t1 and t2
//'
//' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector p_equiv(int m,
                            Rcpp::NumericVector t1, Rcpp::NumericVector t2) {
  if (m < 3) {
    ::Rf_error("Both m must be 3 or greater");
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
  AcceptanceVangel an = AcceptanceVangel(m);
  Rcpp::NumericVector result = Rcpp::NumericVector(num_vals);
  
  for(int i = 0; i < num_vals; ++i) {
    result[i] = an.calc_p_value(t1[i], t2[i]);
  }
  return result;
}



AcceptanceTwoSample::AcceptanceTwoSample(const double n, const double m) :
  AcceptanceBase(m) {
  this->n = n;
  this->dfv_int.init(
      [this](const double v) { return dfv(v); },
      false
  );
  this->dfw_int.init(
      [this](const double w) { return dfw(w); },
      +1, 0., false // integration from 0 to +Inf without oversampling
  );
}

void AcceptanceTwoSample::calculate_factors(const double alpha) {
  auto f = [this, alpha](const double r1) {
    const double cpi_val = cpi(r1);
    const double cpm_val = cpi_val;
    const double r2 = calc_r2(cpm_val);
    const double fjoint = calc_f_joint(r1, r2);
    return cpi_val + cpm_val - fjoint - alpha;
  };
  
  bisection(f, 2, 5, &k1, 500);
  const double cpi_val = cpi(k1);
  k2 = calc_r2(cpi_val);
}

double AcceptanceTwoSample::calc_p_value(const double r1, const double r2) {
  const double cpi_val = cpi(r1);
  const double cpm_val = R::pt(
    r2 / sqrt(1. / this->n + 1. / this->m),
    this->n - 1,
    0,  // lower.tail = false
    0   // log.p = false
  );
  const double cpjoint = calc_f_joint(r1, r2);
  return cpi_val + cpm_val - cpjoint;
}

double AcceptanceTwoSample::dfw(const double w) {
  const double k = n - 1;
  
  return exp(
    (k / 2) * log(k)
    + (k - 1) * log(w)
    - k * pow(w, 2) / 2
    - (R::lgammafn(k / 2) + (k / 2 - 1) * log(2.))
  );
}

double AcceptanceTwoSample::dfv(const double v) {
  const double sqn = sqrt(n);
  return DNORM(v * sqn, false) * sqn;
}

double AcceptanceTwoSample::cpi(const double r1) {
  IntegrationMultDblInf outer_int = IntegrationMultDblInf(
    [this](const double v) { return dfv(v); },
    [r1, this](const double v) {
      IntegrationMultOneInf inner_int = IntegrationMultOneInf(
        [this](const double w) {
          return dfw(w);
        },
        [r1, v, this](const double w) {
          return (1. - pow(
              PNORM(v - r1 * w, false, false),
              m)
          );
        },
        &dfw_int,
        +1, 0.
      );
      return inner_int.result;
    },
    &dfv_int
  );
  return outer_int.result;
}

double AcceptanceTwoSample::calc_r2(const double cpi_val) {
  const double result = R::qt(cpi_val, n - 1., false, false)
    * sqrt(1. / m + 1. / n);
  return result;
}

double AcceptanceTwoSample::calc_f_joint(const double r1, const double r2) {
  IntegrationMultDblInf outer_int = IntegrationMultDblInf(
    [this](const double v) { return dfv(v); },
    [r1, r2, this](const double v) {
      IntegrationMultOneInf inner_int = IntegrationMultOneInf(
        [this](const double w) {
          return dfw(w);
        },
        [r1, r2, v, this](const double w) {
          return calc_f_joint_vangel(v - r1 * w, v - r2 * w);
        },
        &dfw_int,
        +1, 0.
      );
      return inner_int.result;
    },
    &dfv_int
  );
  return outer_int.result;
}

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

TEST_CASE("AcceptanceSample") {
  SUBCASE("dfw & dfv, n=10") {
    AcceptanceTwoSample an = AcceptanceTwoSample(10, 5);
    CHECK_ALMOST_EQ(an.dfw(0.5), 0.1896797, 1e-6);
    CHECK_ALMOST_EQ(an.dfw(1), 1.661563, 1e-6);
    CHECK_ALMOST_EQ(an.dfw(2), 0.0005831514, 1e-6);
    
    CHECK_ALMOST_EQ(an.dfv(0), 1.261566, 1e-6);
    CHECK_ALMOST_EQ(an.dfv(0.5), 0.3614448, 1e-6);
    CHECK_ALMOST_EQ(an.dfv(-0.5), 0.3614448, 1e-6);
    CHECK_ALMOST_EQ(an.dfv(1), 0.008500367, 1e-6);
    CHECK_ALMOST_EQ(an.dfv(-1), 0.008500367, 1e-6);
    CHECK_ALMOST_EQ(an.dfv(2), 2.600282e-09, 1e-9);
    CHECK_ALMOST_EQ(an.dfv(-2), 2.600282e-09, 1e-9);
  }
  SUBCASE("dfw & dfv, n=20") {
    AcceptanceTwoSample an = AcceptanceTwoSample(20, 5);
    CHECK_ALMOST_EQ(an.dfw(0.5), 0.01155585, 1e-6);
    CHECK_ALMOST_EQ(an.dfw(1), 2.437775, 1e-6);
    CHECK_ALMOST_EQ(an.dfw(2), 2.680037e-07, 1e-10);
    
    CHECK_ALMOST_EQ(an.dfv(0.5), 0.1464498, 1e-6);
    CHECK_ALMOST_EQ(an.dfv(-0.5), 0.1464498, 1e-6);
    CHECK_ALMOST_EQ(an.dfv(1), 8.099911e-05, 1e-10);
  }
  SUBCASE("cpi, n=18, m=5") {
    AcceptanceTwoSample an = AcceptanceTwoSample(18, 5);
    CHECK_ALMOST_EQ(an.cpi(2.605), 0.05008773, 1e-6);
  }
  SUBCASE("cpi, n=5, m=18") {
    AcceptanceTwoSample an = AcceptanceTwoSample(5, 18);
    CHECK_ALMOST_EQ(an.cpi(2.605), 0.2946645, 1e-6);
  }
  SUBCASE("factors match prototype R code") {
    AcceptanceTwoSample an = AcceptanceTwoSample(18, 5);
    an.calculate_factors(0.05);
    CHECK_ALMOST_EQ(an.k1, 2.867903, 1e-3);
    CHECK_ALMOST_EQ(an.k2, 1.019985, 1e-3);
  }
  SUBCASE("p-value matches prototype R code") {
    AcceptanceTwoSample an = AcceptanceTwoSample(18, 5);
    const double p = an.calc_p_value(2.867903, 1.019985);
    CHECK_ALMOST_EQ(p, 0.05, 1e-6);
  }
}
