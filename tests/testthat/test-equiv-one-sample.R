suppressMessages(library(dplyr))

test_that("k-factor warnings and errors are raised", {
  expect_error(k_equiv_one_sample(-0.01, 3))
  expect_error(k_equiv_one_sample(1.01, 3))
  expect_warning(k_equiv_one_sample(1e-9, 3))
  expect_warning(k_equiv_one_sample(0.75, 3))
  expect_error(k_equiv_one_sample(0.05, 1))
})

test_that("k-factors match those published in literature", {

  if (requireNamespace("tidyr", quietly = TRUE)) {

    # the files k1.vangel.csv and k2.vangel.csv contain the k factors
    # published in Vangel's 2002 paper.

    k1_vangel <- read.csv(system.file("extdata", "k1.vangel.csv",
                                      package = "cmstatr")) %>%
      tidyr::gather(n, k1, X2:X10) %>%
      mutate(n = as.numeric(substring(n, 2)))

    k2_vangel <- read.csv(system.file("extdata", "k2.vangel.csv",
                                      package = "cmstatr")) %>%
      tidyr::gather(n, k2, X2:X10) %>%
      mutate(n = as.numeric(substring(n, 2)))

    diff_df <- inner_join(k1_vangel, k2_vangel, by = c("alpha", "n")) %>%
      mutate(error_threshold = case_when(n == 2 & alpha > 0.25 ~ 0.10,
                                         n == 2 & alpha <= 0.25 ~ 0.02,
                                         n > 2 & alpha > 0.25 ~ 0.005,
                                         TRUE ~ 1e-3))

    diff_df <- diff_df %>%
      tidyr::gather(fct, vangel, k1:k2) %>%
      group_by(alpha, n) %>%
      filter(n >= 3) %>%
      mutate(calc = k_equiv_one_sample(first(alpha), first(n))) %>%
      mutate(diff = (vangel - calc) / vangel) %>%
      ungroup() %>%
      rowwise() %>%
      mutate(check = expect_lte(abs(diff), error_threshold, label =
                                  paste0("Validation failure for ",
                                         "alpha=", alpha,
                                         ", n=", n,
                                         " computed ", fct, "=", calc,
                                         " but validation ", fct, "=", vangel,
                                         ".\n")))
  }
})
