#include "doctest.h"
#include <cmath>

#define CHECK_ALMOST_EQ(a, b, tol) CHECK(fabs(a - b) <= tol)
