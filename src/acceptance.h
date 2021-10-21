#ifndef C_CODE_ACCEPTANCE_H
#define C_CODE_ACCEPTANCE_H

#include <functional>
#include "integration.h"

#define PNORM(x,lt,lg) R::pnorm5(x,0,1,lt,lg)
#define DNORM(x,lg) R::dnorm4(x,0,1,lg)
#define QNORM(p,lt,lg) R::qnorm5(p,0,1,lt,lg)

class AcceptanceBase {
public:
  AcceptanceBase(const double m, const double alpha);
  double calc_f_joint_vangel(const double t1, const double t2);
  
protected:
  double calc_lambda(const double t1,
                     const double t2, const double x0);
  double h(const double t);
  double a_fcn(const double t);
  double m;
  double alpha;
  IntegrationDblInf a_int;
};

class AcceptanceVangel :
  public AcceptanceBase {
public:
  AcceptanceVangel(const double m, const double alpha);
  double calc_f_min(const double t1);
  double calc_f_mean(const double t2);
  
public:
  double k1;
  double k2;

};

class AcceptanceNew :
  public AcceptanceBase {
  
public:
  AcceptanceNew(const double n, const double m,
                const double alpha, const bool skip_computation = false);
  
  double dfw(const double w);
  double dfv(const double v);
  double cpi(const double r1);
  double calc_r2(const double cpi_val);
  double calc_f_joint(const double r1, const double r2);
  
public:
  double k1;
  double k2;
  
protected:
  double n;
};

#endif