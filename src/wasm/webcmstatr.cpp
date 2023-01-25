//#define CATCH_AMALGAMATED_CUSTOM_MAIN 1
//#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include <stdio.h>
#include <emscripten/emscripten.h>


unsigned int Factorial( unsigned int number ) {
  return number <= 1 ? number : Factorial(number-1)*number;
}

TEST_CASE( "Factorials are computed", "[factorial]" ) {
  REQUIRE( Factorial(1) == 1 );
  REQUIRE( Factorial(2) == 2 );
  REQUIRE( Factorial(3) == 6 );
  REQUIRE( Factorial(10) == 3628800 );
}


#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN
#endif
EXTERN EMSCRIPTEN_KEEPALIVE int do_tests() {
  Catch::Session session; // There must be exactly one instance
  int numFailed = session.run();
  return numFailed;
}

int main() {
  /*const auto f = [](const double x) {
   return cos(x);
   };
   Integration res = Integration(f, 0, M_PI / 2.);
   printf("Integral = %f, should be 1.0\n", res.result);*/
#ifndef WASM
  printf("Hello World\n");
#else
  printf("Hello WASM!\n");
#endif
  /*
  printf("qnorm(.99999) = %f\n", nmath::qnorm5(-0.000100005, 0, 1, true, true));
  printf("dnorm(1) = %f\n", nmath::dnorm4(1, 0, 1, false));
  printf("pnorm(0) = %f\n", nmath::pnorm5(0, 0, 1, true, false));

  printf("pt(0, 10) = %f\n", nmath::pt(0, 10, true, false));
  printf("qt(0.5, 10) = %f\n", nmath::qt(0.5, 10, true, false));
  printf("lgammafn(3) = %f\n", nmath::lgammafn(3));

  AcceptanceTwoSample an = AcceptanceTwoSample(18, 5);
  an.calculate_factors(0.05);
  printf("For n=18, m=5, alpha=0.05: k1=%f, k2=%f\n", an.k1, an.k2);
  const double p = an.calc_p_value(2.867903, 1.019985);
  printf("For k1, k2 = 2.867903, 1.019985, two-sample p = %f (s/b 0.05)\n", p);

  AcceptanceVangel av = AcceptanceVangel(5);
  av.calculate_factors(0.05);
  printf("Vangel, m=5, alpha=0.05: k1=%f, k2=%f\n", av.k1, av.k2);
   */
  return 0;
}

