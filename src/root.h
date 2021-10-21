#ifndef C_CODE_ROOT_H
#define C_CODE_ROOT_H

#include <functional>


#define ROOT_RESULT_SUCCESS 0
#define ROOT_RESULT_MAX_ITT 1
#define ROOT_RESULT_X_CHANGE_TOO_SMALL 2
#define ROOT_FAILED_TO_BRACKET_ROOT 3


int root(std::function<double(const double)> const& f,
         std::function<double(const double)> const& f_prime,
         double x0,
         double* root, int max_itt = 100);
int bisection(std::function<double(const double)> const& f,
              double x1, double x2, double *root,
              int max_itt = 100);

#endif  // C_CODE_ROOT_H
