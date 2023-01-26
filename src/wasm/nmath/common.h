#ifndef COMMON_H
#define COMMON_H

#include <cmath>
#include <math.h>
#include <stdio.h>
#include <float.h>

namespace R {

#define Rboolean bool

#define ISNAN(x) (isnan(x)!=0)
#define R_FINITE(x)    isfinite(x)

#define ML_POSINF	(1.0 / 0.0)
#define ML_NEGINF	((-1.0) / 0.0)
#define ML_NAN		(0.0 / 0.0)

/*
#define ML_POSINF	std::numeric_limits<double>::infinity()
#define ML_NEGINF	-std::numeric_limits<double>::infinity()
#define ML_NAN		std::numeric_limits<double>::quiet_NaN()
*/

/* ----- The following constants and entry points are part of the R API ---- */

/* 30 Decimal-place constants */
/* Computed with bc -l (scale=32; proper round) */

/* SVID & X/Open Constants */
/* Names from Solaris math.h */

#ifndef M_E
#define M_E		2.718281828459045235360287471353	/* e */
#endif

#ifndef M_LOG2E
#define M_LOG2E		1.442695040888963407359924681002	/* log2(e) */
#endif

#ifndef M_LOG10E
#define M_LOG10E	0.434294481903251827651128918917	/* log10(e) */
#endif

#ifndef M_LN2
#define M_LN2		0.693147180559945309417232121458	/* ln(2) */
#endif

#ifndef M_LN10
#define M_LN10		2.302585092994045684017991454684	/* ln(10) */
#endif

#ifndef M_PI
#define M_PI		3.141592653589793238462643383280	/* pi */
#endif

#ifndef M_2PI
#define M_2PI		6.283185307179586476925286766559	/* 2*pi */
#endif

#ifndef M_PI_2
#define M_PI_2		1.570796326794896619231321691640	/* pi/2 */
#endif

#ifndef M_PI_4
#define M_PI_4		0.785398163397448309615660845820	/* pi/4 */
#endif

#ifndef M_1_PI
#define M_1_PI		0.318309886183790671537767526745	/* 1/pi */
#endif

#ifndef M_2_PI
#define M_2_PI		0.636619772367581343075535053490	/* 2/pi */
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI	1.128379167095512573896158903122	/* 2/sqrt(pi) */
#endif

#ifndef M_SQRT2
#define M_SQRT2		1.414213562373095048801688724210	/* sqrt(2) */
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2	0.707106781186547524400844362105	/* 1/sqrt(2) */
#endif

/* R-Specific Constants */

#ifndef M_SQRT_3
#define M_SQRT_3	1.732050807568877293527446341506	/* sqrt(3) */
#endif

#ifndef M_SQRT_32
#define M_SQRT_32	5.656854249492380195206754896838	/* sqrt(32) */
#endif

#ifndef M_LOG10_2
#define M_LOG10_2	0.301029995663981195213738894724	/* log10(2) */
#endif

#ifndef M_SQRT_PI
#define M_SQRT_PI	1.772453850905516027298167483341	/* sqrt(pi) */
#endif

#ifndef M_1_SQRT_2PI
#define M_1_SQRT_2PI	0.398942280401432677939946059934	/* 1/sqrt(2pi) */
#endif

#ifndef M_SQRT_2dPI
#define M_SQRT_2dPI	0.797884560802865355879892119869	/* sqrt(2/pi) */
#endif


#ifndef M_LN_2PI
#define M_LN_2PI	1.837877066409345483560659472811	/* log(2*pi) */
#endif

#ifndef M_LN_SQRT_PI
#define M_LN_SQRT_PI	0.572364942924700087071713675677	/* log(sqrt(pi))
								   == log(pi)/2 */
#endif

#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI	0.918938533204672741780329736406	/* log(sqrt(2*pi))
								 == log(2*pi)/2 */
#endif

#ifndef M_LN_SQRT_PId2
#define M_LN_SQRT_PId2	0.225791352644727432363097614947	/* log(sqrt(pi/2))
								   == log(pi/2)/2 */
#endif



#ifdef HAVE_NEARBYINT
# define R_forceint(x)   nearbyint(x)
#else
# define R_forceint(x)   round(x)
#endif


# define R_nonint(x) 	  (fabs((x) - R_forceint(x)) > 1e-7*fmax2(1., fabs(x)))


#define ME_NONE		0
/*	no error */
#define ME_DOMAIN	1
/*	argument out of domain */
#define ME_RANGE	2
/*	value out of range */
#define ME_NOCONV	4
/*	process did not converge */
#define ME_PRECISION	8
/*	does not have "full" precision */
#define ME_UNDERFLOW	16
/*	and underflow occured (important for IEEE)*/

#define ML_WARNING(x, s) { \
   if(x > ME_DOMAIN) { \
       switch(x) { \
       case ME_DOMAIN: \
	   printf("argument out of domain in '%s'\n", s);	\
	   break; \
       case ME_RANGE: \
	   printf("value out of range in '%s'\n", s);	\
	   break; \
       case ME_NOCONV: \
	   printf("convergence failed in '%s'\n", s);	\
	   break; \
       case ME_PRECISION: \
	   printf("full precision may not have been achieved in '%s'\n", s); \
	   break; \
       case ME_UNDERFLOW: \
	   printf("underflow occurred in '%s'\n", s);	\
	   break; \
       } \
   } \
}
#define ML_WARN_return_NAN { ML_WARNING(ME_DOMAIN, ""); return ML_NAN; }

#define MATHLIB_WARNING(x1, x2) printf(x1,x2)
#define MATHLIB_WARNING2(x1, x2, x3) printf(x1,x2,x3)
#define MATHLIB_WARNING3(x1, x2, x3, x4) printf(x1,x2,x3,x4)
#define MATHLIB_WARNING4(x1, x2, x3, x4, x5) printf(x1,x2,x3,x4,x5)
#define MATHLIB_WARNING5(x1, x2, x3, x4, x5, x6) printf(x1,x2,x3,x4,x5, x6)


int Rf_i1mach(int i);
double Rf_d1mach(int i);
void bratio(double a, double b, double x, double y, double *w, double *w1,
       int *ierr, int log_p);
double fmin2(double x, double y);
double fmax2(double x, double y);
double cospi(double x);
double sinpi(double x);
double tanpi(double x);
double lgammacor(double x);
double chebyshev_eval(double x, const double *a, const int n);
double stirlerr(double n);
double bd0(double x, double np);
void ebd0(double x, double M, double *yh, double *yl);
double dpois_raw(double x, double lambda, int give_log);
double dpois(double x, double lambda, int give_log);
double log1pmx (double x);
double logspace_add (double logx, double logy);

double norm_rand(void);

typedef enum {
    BUGGY_KINDERMAN_RAMAGE,
    AHRENS_DIETER,
    BOX_MULLER,
    USER_NORM,
    INVERSION,
    KINDERMAN_RAMAGE
} N01type;


} // end namespace

#endif  // COMMON_H