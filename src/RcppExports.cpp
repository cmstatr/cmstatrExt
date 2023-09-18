// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// iso_equiv_two_sample
DataFrame iso_equiv_two_sample(const int n, const int m, const double alpha, double t1max, double t2max, const double n_points);
RcppExport SEXP _cmstatrExt_iso_equiv_two_sample(SEXP nSEXP, SEXP mSEXP, SEXP alphaSEXP, SEXP t1maxSEXP, SEXP t2maxSEXP, SEXP n_pointsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type t1max(t1maxSEXP);
    Rcpp::traits::input_parameter< double >::type t2max(t2maxSEXP);
    Rcpp::traits::input_parameter< const double >::type n_points(n_pointsSEXP);
    rcpp_result_gen = Rcpp::wrap(iso_equiv_two_sample(n, m, alpha, t1max, t2max, n_points));
    return rcpp_result_gen;
END_RCPP
}
// k_equiv_two_sample
Rcpp::NumericVector k_equiv_two_sample(double alpha, int n, int m);
RcppExport SEXP _cmstatrExt_k_equiv_two_sample(SEXP alphaSEXP, SEXP nSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(k_equiv_two_sample(alpha, n, m));
    return rcpp_result_gen;
END_RCPP
}
// p_equiv
Rcpp::NumericVector p_equiv(int m, Rcpp::NumericVector t1, Rcpp::NumericVector t2);
RcppExport SEXP _cmstatrExt_p_equiv(SEXP mSEXP, SEXP t1SEXP, SEXP t2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type t1(t1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type t2(t2SEXP);
    rcpp_result_gen = Rcpp::wrap(p_equiv(m, t1, t2));
    return rcpp_result_gen;
END_RCPP
}
// p_equiv_two_sample
Rcpp::NumericVector p_equiv_two_sample(int n, int m, Rcpp::NumericVector t1, Rcpp::NumericVector t2);
RcppExport SEXP _cmstatrExt_p_equiv_two_sample(SEXP nSEXP, SEXP mSEXP, SEXP t1SEXP, SEXP t2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type t1(t1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type t2(t2SEXP);
    rcpp_result_gen = Rcpp::wrap(p_equiv_two_sample(n, m, t1, t2));
    return rcpp_result_gen;
END_RCPP
}
// power_sim_dual_generic
DataFrame power_sim_dual_generic(const int n_qual, const int m_equiv, NumericVector replicates, std::string distribution, Function dist_function, DataFrame param_qual, DataFrame param_equiv, const double k1, const double k2);
RcppExport SEXP _cmstatrExt_power_sim_dual_generic(SEXP n_qualSEXP, SEXP m_equivSEXP, SEXP replicatesSEXP, SEXP distributionSEXP, SEXP dist_functionSEXP, SEXP param_qualSEXP, SEXP param_equivSEXP, SEXP k1SEXP, SEXP k2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n_qual(n_qualSEXP);
    Rcpp::traits::input_parameter< const int >::type m_equiv(m_equivSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type replicates(replicatesSEXP);
    Rcpp::traits::input_parameter< std::string >::type distribution(distributionSEXP);
    Rcpp::traits::input_parameter< Function >::type dist_function(dist_functionSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type param_qual(param_qualSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type param_equiv(param_equivSEXP);
    Rcpp::traits::input_parameter< const double >::type k1(k1SEXP);
    Rcpp::traits::input_parameter< const double >::type k2(k2SEXP);
    rcpp_result_gen = Rcpp::wrap(power_sim_dual_generic(n_qual, m_equiv, replicates, distribution, dist_function, param_qual, param_equiv, k1, k2));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP run_testthat_tests(void *);

static const R_CallMethodDef CallEntries[] = {
    {"_cmstatrExt_iso_equiv_two_sample", (DL_FUNC) &_cmstatrExt_iso_equiv_two_sample, 6},
    {"_cmstatrExt_k_equiv_two_sample", (DL_FUNC) &_cmstatrExt_k_equiv_two_sample, 3},
    {"_cmstatrExt_p_equiv", (DL_FUNC) &_cmstatrExt_p_equiv, 3},
    {"_cmstatrExt_p_equiv_two_sample", (DL_FUNC) &_cmstatrExt_p_equiv_two_sample, 4},
    {"_cmstatrExt_power_sim_dual_generic", (DL_FUNC) &_cmstatrExt_power_sim_dual_generic, 9},
    {"run_testthat_tests", (DL_FUNC) &run_testthat_tests, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_cmstatrExt(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
