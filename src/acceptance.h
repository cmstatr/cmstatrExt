#ifndef C_CODE_ACCEPTANCE_H
#define C_CODE_ACCEPTANCE_H

#include <functional>
#include "integration.h"

#define PNORM(x,lt,lg) R::pnorm5(x,0,1,lt,lg)
#define DNORM(x,lg) R::dnorm4(x,0,1,lg)
#define QNORM(p,lt,lg) R::qnorm5(p,0,1,lt,lg)

class AcceptanceBase {
public:
  AcceptanceBase(const double m);
  double calc_f_joint_vangel(const double t1, const double t2) const;
  
protected:
  double calc_lambda(const double t1,
                     const double t2, const double x0) const;
  double h(const double t) const;
  double a_fcn(const double t) const;
  double m;
  IntegrationDblInf a_int;
};

class AcceptanceVangel :
  public AcceptanceBase {
public:
  AcceptanceVangel(const double m);
  double calc_f_min(const double t1) const;
  double calc_f_mean(const double t2) const;
  void calculate_factors(const double alpha);
  double calc_p_value(const double r1, const double r2) const;
  
public:
  double k1;
  double k2;

};

class AcceptanceTwoSample :
  public AcceptanceBase {
  
public:
  AcceptanceTwoSample(const double n, const double m);
  
  double dfw(const double w) const;
  double dfv(const double v) const;
  double cpi(const double r1) const;
  double calc_r1(const double cpi_val) const;
  double calc_r2(const double cpm_val) const;
  double calc_f_joint(const double r1, const double r2) const;
  void calculate_factors(const double alpha);
  double calc_p_value(const double r1, const double r2) const;
  
public:
  double k1;
  double k2;
  
protected:
  double n;
  IntegrationDblInf dfv_int{};
  IntegrationOneInf dfw_int{};
};

#endif