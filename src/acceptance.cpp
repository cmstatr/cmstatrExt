#include <Rcpp.h>
#include <cmath>
#include "root.h"
#include "integration.h"
#include "acceptance.h"
#include "doctest.h"
#include "doctest-ext.h"


AcceptanceBase::AcceptanceBase(const double m, const double alpha) {
  this->m = m;
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
    return (m - 1.) / m * (h(lam) - lam) + t2 - t1;
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
    const int retval_bisection = bisection(f, -1000, 1000, &result);
    if (retval_bisection != ROOT_RESULT_SUCCESS) {
      ::Rf_error("Root failed. (Newton code=%i, bisection code=%i)",
                 retval, retval_bisection);
    }
  }
  
  return result;
}

double AcceptanceVangel::calc_f_joint(const double t1, const double t2) {
  const auto a_m = [this](const double t) { return a_fcn(t); };
  const auto f1 = [this, t2](const double t) {
    return PNORM(-sqrt(m) * t2, true, false);
  };
  const auto f2 = [this, t1](const double t) {
    return PNORM(
      sqrt(m) * (-t1 + (m - 1.) / m * (h(t) - t)),
      true, false
    );
  };
  
  const double lam = calc_lambda(t1, t2, 0.);  // TODO: Change x0
  
  const IntegrationMultInf num1 = IntegrationMultInf(a_m, f1, &a_int, -1., lam);
  const IntegrationMultInf num2 = IntegrationMultInf(a_m, f2, &a_int, +1., lam);
  
  return (num1.result + num2.result) / a_int.result;
}

double AcceptanceVangel::calc_f_min(const double t1) {
  return 1. - pow(PNORM(-t1, false, false), m);
}

double AcceptanceVangel::calc_f_mean(const double t2) {
  return PNORM(-sqrt(m) * t2, true, false);
}

AcceptanceVangel::AcceptanceVangel(const double m, const double alpha)
  : AcceptanceBase(m, alpha),
    a_int(IntegrationDblInf([this](const double t) { return a_fcn(t); }, true)) {
  
  auto calc_t2 = [this](const double t1) {
    return -QNORM(
        1. - pow(PNORM(-t1, false, false), this->m),
        true, false
    ) / sqrt(this->m);
  };
  
  auto f = [this, calc_t2, alpha](const double t1) {
    double t2 = calc_t2(t1);
    const double fx1 = calc_f_min(t1);
    const double fxbar = calc_f_mean(t2);
    const double fjoint = calc_f_joint(t1, t2);
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

/*
AcceptanceNew::AcceptanceNew(const double n, const double m,
                             const double alpha) :
  AcceptanceBase(m, alpha) {
  this->n = n;
}

double AcceptanceNew::dfw(const double w) {
  const double k = n - 1;
  
  return pow(k, k / 2) * pow(w, k - 1) *
    exp(-k * pow(w, 2) / 2) / (
        gamma(k / 2) * pow(2., k / 2 - 1)
    );
}

double AcceptanceNew::dfv(const double v) {
  return dnorm(v * sqrt(n), false) * sqrt(n);
}

double AcceptanceNew::cpi(const double r1) {
  IntegrationDblInf outer_int = IntegrationDblInf(
    [r1, this](const double v) {
      IntegrationOneInf inner_int = IntegrationOneInf(
        [r1, v, this](const double w) {
          return (1. - pow(
              pnorm(v - r1 * w, false, false),
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
*/