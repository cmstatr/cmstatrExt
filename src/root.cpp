#include "root.h"
#include <cmath>
#include <cfloat>
#include "doctest.h"
#include "doctest-ext.h"


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

TEST_CASE("root") {
  double result = 0.;
  auto f = [](const double x) { return x * x - 4; };
  auto f_prime = [](const double x) { return 2. * x; };
  
  root(f, f_prime, 0.5, &result);
  CHECK_ALMOST_EQ(result, 2., 1e-6);
  
  root(f, f_prime, -0.5, &result);
  CHECK_ALMOST_EQ(result, -2., 1e-6);
}

int bisection(std::function<double(const double)> const& f,
              double x1, double x2, double* root, int max_itt) {
  
  const double abstol = pow(DBL_EPSILON, 0.25);
  int i;
  double f1 = f(x1);
  double f2 = f(x2);
  double fm;
  
  if (f1 * f2 > 0) {
    throw("Endpoints do not have opposite signs");
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

TEST_CASE("bisection") {
  double result;
  auto f = [](const double x) { return x * x - 4; };
  
  bisection(f, 0.5, 10., &result, 100);
  CHECK_ALMOST_EQ(result, 2., 1e-5);
  
  bisection(f, -0.5, -10., &result, 100);
  CHECK_ALMOST_EQ(result, -2., 1e-5);
}
