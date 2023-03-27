test_that("equiv_sample", {
  k <- k_equiv_two_sample(0.05, 18, 5)
  expect_equal(k[1], 2.867903, tolerance = 1e-3)
  expect_equal(k[2], 1.019985, tolerance = 1e-3)

  # the following used to fail
  k_equiv_two_sample(0.05, 100, 6)
})

test_that("equiv_sample simulation", {
  skip_on_cran()

  k <- k_equiv_two_sample(0.05, 30, 9)
  set.seed(1008)
  power_null <- power_sim_dual(30, 9, 5000, "rnorm",
                               data.frame(mean = 0, sd = 1),
                               data.frame(mean = 0, sd = 1),
                               k[1], k[2])
  expect_lt(abs(power_null[["Rejection Rate"]] - 0.05), 0.003)

  k[1] <- k[1] + 1
  p <- p_equiv_two_sample(30, 9, k[1], k[2])
  set.seed(1009)
  power_null <- power_sim_dual(30, 9, 5000, "rnorm",
                               data.frame(mean = 0, sd = 1),
                               data.frame(mean = 0, sd = 1),
                               k[1], k[2])
  expect_lt(abs(power_null[["Rejection Rate"]] - p), 0.003)

  k[1] <- k[1] - 1
  k[2] <- k[2] - 0.3
  p <- p_equiv_two_sample(30, 9, k[1], k[2])
  set.seed(1011)
  power_null <- power_sim_dual(30, 9, 5000, "rnorm",
                               data.frame(mean = 0, sd = 1),
                               data.frame(mean = 0, sd = 1),
                               k[1], k[2])
  expect_lt(abs(power_null[["Rejection Rate"]] - p), 0.003)
})
