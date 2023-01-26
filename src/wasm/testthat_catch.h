// This file contains macros to allow testthat tests to work with catch2

#define context(x)      TEST_CASE(x)
#define test_that(x)    SECTION(x)
#define expect_true(x)  REQUIRE(x)
