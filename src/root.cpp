#include "root.h"
#include <cmath>
#include <cfloat>

#ifndef WASM
#include <testthat.h>
#else // WASM
#include "wasm/catch.hpp"
#include "wasm/testthat_catch.h"
#endif // WASM

#include "testthat-exp.h"


int root(std::function<double(const double)> const& f,
         std::function<double(const double)> const& f_prime,
         double x0,
         double* root, int max_itt) {
  
  const double abstol = pow(DBL_EPSILON, 0.25);
  int i;
  double f0, f_prime_0, x1;
  
  for(i=0; i < max_itt; ++i) {
    f0 = f(x0);
    if (fabs(f0) <= abstol)
    {
      *root = x0;
      return ROOT_RESULT_SUCCESS;
    }
    f_prime_0 = f_prime(x0);
    x1 = x0 - f0 / f_prime_0;
    if (fabs(x1 - x0) <= abstol)
    {
      *root = x1;
      return ROOT_RESULT_X_CHANGE_TOO_SMALL;
    }
    x0 = x1;
  }
  
  return ROOT_RESULT_MAX_ITT;
}

context("root") {
  test_that("root") {
    double result = 0.;
    auto f = [](const double x) { return x * x - 4; };
    auto f_prime = [](const double x) { return 2. * x; };
    
    root(f, f_prime, 0.5, &result);
    expect_almost_eq(result, 2., 1e-6);
    
    root(f, f_prime, -0.5, &result);
    expect_almost_eq(result, -2., 1e-6);
  }
}

int bisection(std::function<double(const double)> const& f,
              double x1, double x2, double* root, int max_itt) {
  
  const double abstol = pow(DBL_EPSILON, 0.25);
  int i;
  double f1 = f(x1);
  double f2 = f(x2);
  double fm;
  
  double delta_1 = 0.05 * fmax(1e-4, fabs(x1));
  double delta_2 = 0.05 * fmax(1e-4, fabs(x2));
  
  for(i = 0; f1 * f2 > 0 && i < max_itt; ++i) {
    x1 -= delta_1;
    x2 += delta_2;
    f1 = f(x1);
    f2 = f(x2);
    delta_1 *= 2;
    delta_2 *= 2;
  }
  if (i == max_itt) {
    return ROOT_FAILED_TO_BRACKET_ROOT;
  }
  if (fabs(f1) <= abstol) {
    *root = x1;
    return ROOT_RESULT_SUCCESS;
  }
  if (fabs(f2) <= abstol) {
    *root = x2;
    return ROOT_RESULT_SUCCESS;
  }
  
  for(i = 0; i < max_itt; ++i) {
    *root = (x1 + x2) / 2.;
    fm = f(*root);
    if (fabs(fm) <= abstol) {
      return ROOT_RESULT_SUCCESS;
    }
    if (f1 * fm < 0) {
      x2 = *root;
      f2 = fm;
    } else {
      x1 = *root;
      f1 = fm;
    }
  }
  return ROOT_RESULT_MAX_ITT;
}

context("bisection") {
  test_that("bisection") {
    double result;
    auto f = [](const double x) { return x * x - 4; };
    
    bisection(f, 0.5, 10., &result, 100);
    expect_almost_eq(result, 2., 1e-5);
    
    bisection(f, -0.5, -10., &result, 100);
    expect_almost_eq(result, -2., 1e-5);
  }
}
