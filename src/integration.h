#ifndef C_CODE_INTEGRATION_H
#define C_CODE_INTEGRATION_H

#include <functional>

#define INT_RESULT_SUCCESS 0
#define INT_RESULT_SUB_TOO_SMALL 1
#define INT_RESULT_MAX_SEGMENTS 2
#define NK 15
#define MAX_SEGMENTS 100

class Segment {
public:
  double result{};
  double resabs{};
  double resasc{};
  double abserr{};
  double a{};
  double b{};
  double x[NK]{};
  double fv[NK]{};
};

class IntegrationBase {
public:
  int message{};
  double result{};
  double abserr{};
  int num_segments{};
  Segment segments[MAX_SEGMENTS]{};
  
protected:
  static constexpr const double xgk[] = {
    -0.991455371120812639206854697526329,
    -0.949107912342758524526189684047851,
    -0.864864423359769072789712788640926,
    -0.741531185599394439863864773280788,
    -0.586087235467691130294144838258730,
    -0.405845151377397166906606412076961,
    -0.207784955007898467600689403773245,
    0.000000000000000000000000000000000,
    0.207784955007898467600689403773245,
    0.405845151377397166906606412076961,
    0.586087235467691130294144838258730,
    0.741531185599394439863864773280788,
    0.864864423359769072789712788640926,
    0.949107912342758524526189684047851,
    0.991455371120812639206854697526329
  };
  constexpr static const double wgk[] = {
    0.022935322010529224963732008058970,
    0.063092092629978553290700663189204,
    0.104790010322250183839876322541518,
    0.140653259715525918745189590510238,
    0.169004726639267902826583426598550,
    0.190350578064785409913256402421014,
    0.204432940075298892414161999234649,
    0.209482141084727828012999174891714,
    0.204432940075298892414161999234649,
    0.190350578064785409913256402421014,
    0.169004726639267902826583426598550,
    0.140653259715525918745189590510238,
    0.104790010322250183839876322541518,
    0.063092092629978553290700663189204,
    0.022935322010529224963732008058970
  };
  constexpr static const double wg[] = {
    0,
    0.129484966168869693270611432679082,
    0,
    0.279705391489276667901467771423780,
    0,
    0.381830050505118944950369775488975,
    0,
    0.417959183673469387755102040816327,
    0,
    0.381830050505118944950369775488975,
    0,
    0.279705391489276667901467771423780,
    0,
    0.129484966168869693270611432679082,
    0
  };
  
  double total_error() const;
  double total_area() const;
  int max_abserr_segment() const;
  static double rescale_error(double err, double result_abs,
                              double result_asc);
  static bool subdivision_too_small(double a1, double a2, double b2);
  static void integration_qk(const std::function<double(const double)>& f,
                             double a_seg, double b_seg, Segment *res);
  void adapt_quadrature(std::function<double(const double)> const& f);
  void oversample_quadrature(const std::function<double(const double)> &f);
  void qags(std::function<double(const double)> const& f, bool oversample);
  static void integration_qk_mult(const std::function<double(const double)>& g,
                                  const Segment *orig_seg, Segment *new_seg);
  void qags_mult(const std::function<double(const double)>& f,
                 const std::function<double(const double)>& g,
                 const double a, const double b,
                 const IntegrationBase *f_result);
};

class Integration : public IntegrationBase {
public:
  Integration(std::function<double(const double)> const& f, double a, double b,
              bool oversample = false);
  
};

class IntegrationOneInf : public IntegrationBase {
public:
  void init(std::function<double(const double)> const& f, int inf_side,
            double c, bool oversample = false);
};

class IntegrationDblInf : public IntegrationBase {
public:
  void init(std::function<double(const double)> const& f,
            bool oversample = false);
};

class IntegrationMult : public IntegrationBase {
public:
  IntegrationMult(const std::function<double(const double)>& f,
                  const std::function<double(const double)>& g,
                  const Integration *f_result, const double a, const double b);
};

class IntegrationMultOneInf : public IntegrationBase {
public:
  IntegrationMultOneInf(const std::function<double(const double)>& f,
                        const std::function<double(const double)>& g,
                        const IntegrationBase *f_result, const int inf_side,
                        const double c);
};

class IntegrationMultDblInf : public IntegrationBase {
public:
  IntegrationMultDblInf(const std::function<double(const double)>& f,
                        const std::function<double(const double)>& g,
                        const IntegrationBase *f_result);
};

#endif //C_CODE_INTEGRATION_H
