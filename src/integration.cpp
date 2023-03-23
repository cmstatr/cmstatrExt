// [[Rcpp::plugins("cpp17")]]

#include "integration.h"
#include <cmath>
#include <cfloat>

#ifndef WASM
#include <Rcpp.h>
#include <testthat.h>
#else // WASM
#include "wasm/nmath/nmath.h"
#include "wasm/catch.hpp"
#include "wasm/testthat_catch.h"
#endif // WASM

#include "testthat-exp.h"

double IntegrationBase::rescale_error(double err, double result_abs,
                                      double result_asc) {
  double scale;
  double min_err;
  
  err = fabs(err);
  if (result_asc != 0 && err != 0) {
    scale = pow(200 * err / result_asc,  1.5);
    if (scale < 1) {
      err = result_asc * scale;
    } else {
      err = result_asc;
    }
  }
  
  if (result_abs > DBL_MIN / (50 * DBL_EPSILON)) {
    min_err = 50 * DBL_EPSILON * result_abs;
    if (min_err > err) {
      err = min_err;
    }
  }
  
  return err;
}

bool IntegrationBase::subdivision_too_small(double a1, double a2, double b2) {
  const double tmp = (1 + 100 * DBL_EPSILON) * (fabs(a2) + 1000 * DBL_MIN);
  return fabs(a1) <= tmp && fabs(b2) <= tmp;
}

double IntegrationBase::total_error() const {
  double err = 0;
  for(int i = 0; i < num_segments; ++i) {
    err += segments[i].abserr;
  }
  return err;
}

double IntegrationBase::total_area() const {
  double area = 0;
  for(int i = 0; i < num_segments; ++i) {
    area += segments[i].result;
  }
  return area;
}

int IntegrationBase::max_abserr_segment() const {
  double max_err = 0;
  int worst = 0;
  
  for(int i = 0; i < num_segments; ++i) {
    const Segment *s = &segments[i];
    if(s->abserr > max_err) {
      max_err = s->abserr;
      worst = i;
    }
  }
  return worst;
}

void IntegrationBase::integration_qk(
    const std::function<double( const double)>& f,
    double a_seg, double b_seg, Segment* res) {
  
  const double center = 0.5 * (a_seg + b_seg);
  const double half_length = 0.5 * (b_seg - a_seg);
  const double abs_half_length = fabs(half_length);
  double result_gauss = 0;
  double result_kronrod = 0;
  double result_abs = 0;
  double result_asc = 0;
  double err;
  
  for(int i = 0; i < NK; ++i){
    res->x[i] = center + half_length * xgk[i];
    res->fv[i] = f(res->x[i]);
    result_gauss += wg[i] * res->fv[i];
    result_kronrod += wgk[i] * res->fv[i];
    result_abs += wgk[i] * fabs(res->fv[i]);
  }
  
  result_gauss *= half_length;
  result_kronrod *= half_length;
  result_abs *= half_length;
  
  for(int i = 0; i < NK; ++i) {
    result_asc += wgk[i] * fabs(res->fv[i] - result_kronrod * 0.5);
  }
  result_asc *= abs_half_length;
  
  err = result_kronrod - result_gauss;
  
  res->result = result_kronrod;
  res->resabs = result_abs;
  res->resasc = result_asc;
  res->abserr = rescale_error(err, result_abs, result_asc);
  res->a = a_seg;
  res->b = b_seg;
}

void IntegrationBase::adapt_quadrature(
    std::function<double(const double)> const& f) {
  
  const double reltol = pow(DBL_EPSILON, 0.25);
  const double abstol = pow(DBL_EPSILON, 0.25);
  double a1, a2, b1, b2;
  int worst_seg_i, i_sub;
  Segment* worst_seg;
  
  for(i_sub = num_segments; i_sub < MAX_SEGMENTS; ++i_sub) {
    if (total_error() <= fmax(abstol, reltol * total_area())) {
      break;
    }
    
    worst_seg_i = max_abserr_segment();
    worst_seg = &segments[worst_seg_i];
    a1 = worst_seg->a;
    b1 = 0.5 * (worst_seg->a + worst_seg->b);
    a2 = b1;
    b2 = worst_seg->b;
    
    if (subdivision_too_small(a1, a2, b2)) {
      message = INT_RESULT_SUB_TOO_SMALL;
      return;
    }
    
    integration_qk(f, a1, b1, &segments[worst_seg_i]);
    integration_qk(f, a2, b2, &segments[num_segments++]);
  }
  if (i_sub >= MAX_SEGMENTS) {
    message = INT_RESULT_MAX_SEGMENTS;
  }
}

void IntegrationBase::oversample_quadrature(
    std::function<double(const double)> const& f) {
  
  double a1, a2, b1, b2;
  int i_sub;
  Segment* s;
  const int initial_num_segments = num_segments;
  
  for(i_sub = 0;
      i_sub < initial_num_segments && i_sub < MAX_SEGMENTS;
      ++i_sub) {
    s = &segments[i_sub];
    a1 = s->a;
    b1 = 0.5 * (s->a + s->b);
    a2 = b1;
    b2 = s->b;
    
    if (!subdivision_too_small(a1, a2, b2)) {
      integration_qk(f, a1, b1, &segments[i_sub]);
      integration_qk(f, a2, b2, &segments[num_segments++]);
    }
  }
  
  if (i_sub >= MAX_SEGMENTS) {
    message = INT_RESULT_MAX_SEGMENTS;
  }
}

void IntegrationBase::qags(
    std::function<double(const double)> const& f, bool oversample) {
  message = INT_RESULT_SUCCESS;
  
  adapt_quadrature(f);
  
  if (oversample) {
    oversample_quadrature(f);
  }
  
  result = total_area();
  abserr = total_error();
}

void IntegrationBase::integration_qk_mult(
    const std::function<double(const double)>& g,
    const Segment* orig_seg, Segment* new_seg) {
  
  const double half_length = 0.5 * (orig_seg->b - orig_seg->a);
  const double abs_half_length = fabs(half_length);
  double gi, result_gauss = 0, result_kronrod = 0, result_abs = 0, result_asc = 0, err;
  
  new_seg->a = orig_seg->a;
  new_seg->b = orig_seg->b;
  
  for (int i = 0; i < NK; ++i) {
    new_seg->x[i] = orig_seg->x[i];
    gi = g(orig_seg->x[i]);
    new_seg->fv[i] = gi * orig_seg->fv[i];
    
    result_gauss += wg[i] * new_seg->fv[i];
    result_kronrod += wgk[i] * new_seg->fv[i];
    result_abs += wgk[i] * fabs(new_seg->fv[i]);
  }
  
  result_gauss *= half_length;
  result_kronrod *= half_length;
  result_abs *= half_length;
  
  for(int i = 0; i < NK; ++i) {
    result_asc += wgk[i] * fabs(new_seg->fv[i] - result_kronrod * 0.5);
  }
  result_asc *= abs_half_length;
  
  err = result_kronrod - result_gauss;
  
  new_seg->result = result_kronrod;
  new_seg->resabs = result_abs;
  new_seg->resasc = result_asc;
  new_seg->abserr = rescale_error(err, result_abs, result_asc);
}

void IntegrationBase::qags_mult(const std::function<double(const double)>& f,
                                const std::function<double(const double)>& g,
                                const double a, const double b,
                                const IntegrationBase *f_result) {
  auto fg = [f, g](double t) { return f(t) * g(t); };
  
  for (int i_orig = 0; i_orig < f_result->num_segments; ++i_orig) {
    const Segment *s_orig = &f_result->segments[i_orig];
    if (a <= s_orig->a && s_orig->b <= b) {
      // segment is fully within new bounds
      integration_qk_mult(g, &f_result->segments[i_orig], &segments[num_segments++]);
    } else if (s_orig->b <= a || s_orig->a >= b) {
      // segment is fully to the left, do nothing
      // OR: segment is fully to the right, do nothing
    } else {
      // segment is partially within the new bound
      integration_qk(
        fg,
        fmax(a, s_orig->a),
        fmin(b, s_orig->b),
        &segments[num_segments++]
      );
      
    }
  }
  
  adapt_quadrature(fg);
  
  result = total_area();
  abserr = total_error();
}

Integration::Integration(std::function<double(const double)> const& f,
                         double a, double b, bool oversample) {
  integration_qk(f, a, b, &segments[0]);
  num_segments = 1;
  
  qags(f, oversample);
}

context("Basic quadrature") {
  auto f_cos = [](const double x) { return cos(x); };
  
  test_that("Basic integration of cos") {
    Integration res = Integration(f_cos, 0, M_PI / 2.);
    expect_almost_eq(res.result, 1, 1e-6);
  }
  
  test_that("Integration of cos with oversampling") {
    Integration res = Integration(f_cos, 0, M_PI / 2., true);
    expect_almost_eq(res.result, 1, 1e-6);
  }
}

void IntegrationOneInf::init(
  std::function<double(const double)> const& f, int inf_side,
  double c, bool oversample) {
  
  auto f_prime = [f](const double t) {
    return f(tan(t)) * pow(cos(t), -2.);
  };
  
  double a, b;
  
  if (inf_side < 0) {
    // From negative infinity
    a = -M_PI / 2.;
    b = atan(c);
  } else {
    // To positive infinity
    a = atan(c);
    b = M_PI / 2.;
  }
  
  integration_qk(f_prime, a, b, &segments[0]);
  num_segments = 1;
  
  qags(f_prime, oversample);
}

context("Quadrature with one infinite bound") {
  auto fn = [](const double t) { return R::dnorm4(t, 0, 1, 0); };
  
  test_that("From -Inf to a negative number") {
    IntegrationOneInf res_m = IntegrationOneInf();
    res_m.init(fn, -1, -0.5);
    expect_almost_eq(res_m.result, 0.3085375, 1e-6);
  }
  test_that("From -Inf to a positive number") {
    IntegrationOneInf res_m = IntegrationOneInf();
    res_m.init(fn, -1, 0.5);
    expect_almost_eq(res_m.result, 0.6914625, 1e-6);
  }
  test_that("From a negative number to Inf") {
    IntegrationOneInf res_m = IntegrationOneInf();
    res_m.init(fn, 1, -0.5);
    expect_almost_eq(res_m.result, 0.6914625, 1e-6);
  }
  test_that("From a positive number to Inf") {
    IntegrationOneInf res_m = IntegrationOneInf();
    res_m.init(fn, 1, 0.5);
    expect_almost_eq(res_m.result, 0.3085375, 1e-6);
  }
}

void IntegrationDblInf::init(
    std::function<double(const double)> const& f, bool oversample) {
  
  auto f_prime = [f](const double t) {
    return f(tan(t)) * pow(cos(t), -2.);
  };
  
  integration_qk(f_prime, -M_PI / 2., M_PI / 2., &segments[0]);
  num_segments = 1;
  
  qags(f_prime, oversample);
}

context("Quadrature with infinite bounds") {
  auto h = [](double t) {
    return exp(R::dnorm4(t, 0, 1, 1) - R::pnorm5(t, 0, 1, 0, 1));
  };

  auto a_fcn = [h](double t, double m) {
    const double hmt = t > 60. ? pow(t, -1) : h(t) - t;
    return exp(
      - (m - 1) * (R::dnorm4(t, 0, 1, 1) - R::pnorm5(t, 0, 1, 0, 1)) +
        pow(m - 1, 2) / (2 * m) * pow(hmt, 2) +
        (m - 1) * t * hmt
    ) * sqrt(t > 60 ? pow(t, -2) : 1 - h(t) * hmt);
  };
  
  auto a_10 = [a_fcn](double t) { return a_fcn(t, 10); };
  
  test_that("Integration of A") {
    IntegrationDblInf res = IntegrationDblInf();
    res.init(a_10);
    expect_almost_eq(res.result, 1038.764, 0.01);
  }
  
  test_that("Integration of A with oversampling") {
    IntegrationDblInf res = IntegrationDblInf();
    res.init(a_10, true);
    expect_almost_eq(res.result, 1038.764, 0.01);
  }
}

IntegrationMult::IntegrationMult(const std::function<double(const double)>& f,
                                 const std::function<double(const double)>& g,
                                 const Integration *f_result,
                                 const double a, const double b) {
  message = f_result->message;
  num_segments = 0;
  
  qags_mult(f, g, a, b, f_result);
}

context("Mult quadrature") {
  auto h = [](double t) {
    return exp(R::dnorm4(t, 0, 1, 1) - R::pnorm5(t, 0, 1, 0, 1));
  };
  
  auto a_fcn = [h](double t, double m) {
    const double hmt = t > 60. ? pow(t, -1) : h(t) - t;
    return exp(
      - (m - 1) * (R::dnorm4(t, 0, 1, 1) - R::pnorm5(t, 0, 1, 0, 1)) +
        pow(m - 1, 2) / (2 * m) * pow(hmt, 2) +
        (m - 1) * t * hmt
    ) * sqrt(t > 60 ? pow(t, -2) : 1 - h(t) * hmt);
  };
  
  auto a_10 = [a_fcn](double t) { return a_fcn(t, 10); };
  Integration res = Integration(a_10, -100, 100);
  
  auto g = [](double t) { return 1.; };
  
  test_that("Basic integration of A") {
    IntegrationMult res2 = IntegrationMult(a_10, g, &res, -100, 100);
    expect_almost_eq(res2.result, 1038.763750, 0.001);
  }
  
  test_that("Partial range") {
    auto sin_lambda = [](double t) { return sin(t); };
    Integration int_sin = Integration(sin_lambda, -M_PI, M_PI, true);
    expect_almost_eq(int_sin.result, 0., 1e-6);  // should be zero
    
    IntegrationMult res2 = IntegrationMult(sin_lambda, g, &int_sin,
                                           -M_PI, M_PI);
    expect_almost_eq(res2.result, 0., 1e-6);  // should be 0 still
    
    res2 = IntegrationMult(sin_lambda, g, &int_sin, 0, M_PI);
    expect_almost_eq(res2.result, 2, 1e-6);  // should be 2
    
    res2 = IntegrationMult(sin_lambda, g, &int_sin, -M_PI, 0);
    expect_almost_eq(res2.result, -2, 1e-6);  // should be -2
    
    res2 = IntegrationMult(sin_lambda, g, &int_sin, 0, 1);
    expect_almost_eq(res2.result, 0.45970, 1e-4);
  }
  
  test_that("Partial range with many segments") {
    auto lambda = [](double t) { return R::pnorm5(t, 0, 1, 1, 0) * t; };
    Integration res1 = Integration(lambda, -100., 100.);
    expect_true(res1.num_segments > 4);
    expect_almost_eq(res1.result, 4999.491, 1);
    
    IntegrationMult res_m = IntegrationMult(lambda, g, &res1, -10, 0);
    expect_almost_eq(res_m.result, -0.25, 1e-6);
  }
}

IntegrationMultOneInf::IntegrationMultOneInf(
  const std::function<double(const double)>& f,
  const std::function<double(const double)>& g,
  const IntegrationBase *f_result, const int inf_side, const double c) {
  
  auto f_prime = [f](const double t) {
    return f(tan(t)) * pow(cos(t), -2.);
  };
  // f_prime includes the sec^2(t) factor, but g_prime does not so that we can
  // multiply these together
  auto g_prime = [g](const double t) {
    return g(tan(t));
  };
  
  message = f_result->message;
  num_segments = 0;
  
  if (inf_side < 0)
  {
    // From negative infinity
    qags_mult(f_prime, g_prime, -M_PI / 2., atan(c), f_result);
  } else {
    // To positive infinity
    qags_mult(f_prime, g_prime, atan(c), M_PI / 2., f_result);
  }
  
}

context("Mult Quadrature Infinte Range") {
  auto h = [](double t) {
    return exp(R::dnorm4(t, 0, 1, 1) - R::pnorm5(t, 0, 1, 0, 1));
  };
  
  auto a_fcn = [h](double t, double m) {
    const double hmt = t > 60. ? pow(t, -1) : h(t) - t;
    return exp(
      - (m - 1) * (R::dnorm4(t, 0, 1, 1) - R::pnorm5(t, 0, 1, 0, 1)) +
        pow(m - 1, 2) / (2 * m) * pow(hmt, 2) +
        (m - 1) * t * hmt
    ) * sqrt(t > 60 ? pow(t, -2) : 1 - h(t) * hmt);
  };
  
  auto a_10 = [a_fcn](double t) { return a_fcn(t, 10); };
  IntegrationDblInf res = IntegrationDblInf();
  res.init(a_10);
  auto g = [](double t) { return 1.; };
  
  expect_almost_eq(res.result, 1038.764, 0.02);
  
  test_that("From -Inf to a negative number") {
    IntegrationMultOneInf res_m = IntegrationMultOneInf(a_10, g, &res, -1, -1.);
    expect_almost_eq(res_m.result, 797.074, 0.02);
  }
  test_that("From -Inf to a positive number") {
    IntegrationMultOneInf res_m = IntegrationMultOneInf(a_10, g, &res, -1, 1.);
    expect_almost_eq(res_m.result, 1037.554, 0.02);
  }
  test_that("From a negative number to Inf") {
    IntegrationMultOneInf res_m = IntegrationMultOneInf(a_10, g, &res, 1, -1.);
    expect_almost_eq(res_m.result, 241.6897, 0.02);
  }
  test_that("From a positive number to Inf") {
    IntegrationMultOneInf res_m = IntegrationMultOneInf(a_10, g, &res, 1, 1.);
    expect_almost_eq(res_m.result, 1.20975, 0.02);
  }
}

IntegrationMultDblInf::IntegrationMultDblInf(
  const std::function<double(const double)>& f,
  const std::function<double(const double)>& g,
  const IntegrationBase *f_result) {
  
  auto f_prime = [f](const double t) {
    return f(tan(t)) * pow(cos(t), -2.);
  };
  // f_prime includes the sec^2(t) factor, but g_prime does not so that we can
  // multiply these together
  auto g_prime = [g](const double t) {
    return g(tan(t));
  };
  
  message = f_result->message;
  num_segments = 0;
  
  qags_mult(f_prime, g_prime, -M_PI / 2., M_PI / 2., f_result);

}

context("Mult Quadrature Double Infinte Range") {
  auto f = [](const double x){ return R::dnorm4(x, 0, 1, 0); };
  IntegrationDblInf f_int = IntegrationDblInf();
  f_int.init(f);
  
  test_that("From a positive number to Inf") {
    auto g = [](const double x){ return 2.; };
    IntegrationMultDblInf res_m = IntegrationMultDblInf(f, g, &f_int);
    expect_almost_eq(res_m.result, 2, 0.0002);
  }
}
