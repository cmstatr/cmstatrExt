#ifndef NMATH_H
#define NMATH_H

namespace R {

double qnorm5(double p, double mu, double sigma, int lower_tail, int log_p);
double dnorm4(double x, double mu, double sigma, int give_log);
double pnorm5(double x, double mu, double sigma, int lower_tail, int log_p);
double rnorm(double mu, double sigma);
double pt(double x, double n, int lower_tail, int log_p);
double qt(double p, double ndf, int lower_tail, int log_p);
double dt(double x, double n, int give_log);
double lgammafn(double x);
double gammafn(double x);
double lbeta(double a, double b);
double pbeta(double x, double a, double b, int lower_tail, int log_p);


} // end namespace

#endif  // NMATH_H
