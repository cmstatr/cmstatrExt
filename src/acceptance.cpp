#include <Rcpp.h>
#include <cmath>
#include "root.h"
#include "integration.h"
#include "acceptance.h"
#include "doctest.h"
#include "doctest-ext.h"


AcceptanceBase::AcceptanceBase(const double m, const double alpha) :
  m(m),
  a_int(IntegrationDblInf([this](const double t) { return a_fcn(t); }, true)) {
  this->alpha = alpha;
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
  
  const IntegrationMultInf num1 = IntegrationMultInf(a_m, f1, &a_int, -1., lam);
  const IntegrationMultInf num2 = IntegrationMultInf(a_m, f2, &a_int, +1., lam);
  
  return (PNORM(sqrt(m) * t2, true, false) * num1.result + num2.result) /
    a_int.result;
}

double AcceptanceVangel::calc_f_min(const double t1) {
  return 1. - pow(PNORM(t1, false, false), m);
}

double AcceptanceVangel::calc_f_mean(const double t2) {
  return PNORM(sqrt(m) * t2, true, false);
}

AcceptanceVangel::AcceptanceVangel(const double m, const double alpha)
  : AcceptanceBase(m, alpha) {
  
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
  
  int retval = bisection(f, -0.1, 1, &k1, 100);
  k2 = calc_t2(k1);
}

TEST_CASE("Acceptance Vangel") {
  SUBCASE("m=5, alpha=0.05") {
    AcceptanceVangel ag = AcceptanceVangel(5, 0.05);
    CHECK_ALMOST_EQ(ag.k1, 2.5286, 0.001);
    CHECK_ALMOST_EQ(ag.k2, 0.8525, 0.001);
  }
  SUBCASE("m=10, alpha=0.05") {
    AcceptanceVangel ag = AcceptanceVangel(10, 0.05);
    CHECK_ALMOST_EQ(ag.k1, 2.7772, 0.001);
    CHECK_ALMOST_EQ(ag.k2, 0.6089, 0.001);
  }
  SUBCASE("m=5, alpha=0.5") {
    AcceptanceVangel ag = AcceptanceVangel(5, 0.5);
    CHECK_ALMOST_EQ(ag.k1, 1.3498, 0.001);
    CHECK_ALMOST_EQ(ag.k2, 0.1473, 0.001);
  }
  SUBCASE("m=10, alpha=0.5") {
    AcceptanceVangel ag = AcceptanceVangel(10, 0.5);
    CHECK_ALMOST_EQ(ag.k1, 1.7258, 0.001);
    CHECK_ALMOST_EQ(ag.k2, 0.1217, 0.001);
  }
  SUBCASE("m=5, alpha=0.0005") {
    AcceptanceVangel ag = AcceptanceVangel(5, 0.0005);
    CHECK_ALMOST_EQ(ag.k1, 3.8864, 0.05);
    CHECK_ALMOST_EQ(ag.k2, 1.5546, 0.05);
  }
  SUBCASE("m=10, alpha=0.0005") {
    AcceptanceVangel ag = AcceptanceVangel(10, 0.0005);
    CHECK_ALMOST_EQ(ag.k1, 4.0541, 0.05);
    CHECK_ALMOST_EQ(ag.k2, 1.1002, 0.05);
  }
}


AcceptanceTwoSample::AcceptanceTwoSample(const double n, const double m,
                             const double alpha) :
  AcceptanceBase(m, alpha) {
  this->n = n;
}

void AcceptanceTwoSample::calculate_factors() {
  auto f = [this](const double r1) {
    const double cpi_val = cpi(r1);
    const double cpm_val = cpi_val;
    const double r2 = calc_r2(cpm_val);
    const double fjoint = calc_f_joint(r1, r2);
    return cpi_val + cpm_val - fjoint - alpha;
  };
  
  int retval = bisection(f, 2, 5, &k1, 500);
  const double cpi_val = cpi(k1);
  k2 = calc_r2(cpi_val);
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
  IntegrationDblInf outer_int = IntegrationDblInf(
    [r1, this](const double v) {
      IntegrationOneInf inner_int = IntegrationOneInf(
        [r1, v, this](const double w) {
          return (1. - pow(
              PNORM(v - r1 * w, false, false),
              m)
          ) * dfw(w);
        },
        +1, 0., false // integration from 0 to +Inf without oversampling
      );
      
      return inner_int.result * dfv(v);
    },
    false // no oversampling
  );
  return outer_int.result;
}

double AcceptanceTwoSample::calc_r2(const double cpi_val) {
  const double result = R::qt(cpi_val, n - 1., false, false)
    * sqrt(1 / m + 1 / n);
  return result;
}

double AcceptanceTwoSample::calc_f_joint(const double r1, const double r2) {
  IntegrationDblInf outer_int = IntegrationDblInf(
    [r1, r2, this](const double v) {
      IntegrationOneInf inner_int = IntegrationOneInf(
        [r1, r2, v, this](const double w) {
          return calc_f_joint_vangel(v - r1 * w, v - r2 * w) * dfw(w);
        },
        +1, 0., false // integration from 0 to +Inf without oversampling
      );
      return inner_int.result * dfv(v);
    },
    false // no oversampling
  );
  return outer_int.result;
}

//' Calculate the factors for a two-sample acceptance test
//'
//' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector k_equiv_two_sample(int n, int m, double alpha) {
  if (n < 3 || m < 3) {
    ::Rf_error("Both n and m must be 3 or greater");
  }
  
  AcceptanceTwoSample an = AcceptanceTwoSample(n, m, alpha);
  an.calculate_factors();
  
  return Rcpp::NumericVector::create(
    an.k1,
    an.k2
  );
}

TEST_CASE("AcceptanceSample") {
  SUBCASE("dfw & dfv, n=10") {
    AcceptanceTwoSample an = AcceptanceTwoSample(10, 5, 0.05);
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
    AcceptanceTwoSample an = AcceptanceTwoSample(20, 5, 0.05);
    CHECK_ALMOST_EQ(an.dfw(0.5), 0.01155585, 1e-6);
    CHECK_ALMOST_EQ(an.dfw(1), 2.437775, 1e-6);
    CHECK_ALMOST_EQ(an.dfw(2), 2.680037e-07, 1e-10);
    
    CHECK_ALMOST_EQ(an.dfv(0.5), 0.1464498, 1e-6);
    CHECK_ALMOST_EQ(an.dfv(-0.5), 0.1464498, 1e-6);
    CHECK_ALMOST_EQ(an.dfv(1), 8.099911e-05, 1e-10);
  }
  SUBCASE("cpi, n=18, m=5") {
    AcceptanceTwoSample an = AcceptanceTwoSample(18, 5, 0.05);
    CHECK_ALMOST_EQ(an.cpi(2.605), 0.05008773, 1e-6);
  }
  SUBCASE("cpi, n=5, m=18") {
    AcceptanceTwoSample an = AcceptanceTwoSample(5, 18, 0.05);
    CHECK_ALMOST_EQ(an.cpi(2.605), 0.2946645, 1e-6);
  }
  SUBCASE("factors match prototype R code") {
    AcceptanceTwoSample an = AcceptanceTwoSample(18, 5, 0.05);
    an.calculate_factors();
    CHECK_ALMOST_EQ(an.k1, 2.867903, 1e-3);
    CHECK_ALMOST_EQ(an.k2, 1.019985, 1e-3);
  }
}
