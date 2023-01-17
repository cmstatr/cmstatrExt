test_that("power_sim_dual custom distribution matches canned distribution", {

  qual <- data.frame(mean = 0, sd = 1)
  eq <- data.frame(mean = c(-2, -1.5, -1, -.5, 0), sd = 1)

  set.seed(100)
  res_canned <- power_sim_dual(18, 6,
                               c(2500, 2500),
                               "norm",
                               qual, eq,
                               2.9594727, 0.9541395)
  rnorm1 <- rnorm
  set.seed(100)
  res_custom <- power_sim_dual(18, 6,
                               2500,  # should be the same as c(2500, 2500)
                               "rnorm1",
                               qual, eq,
                               2.9594727, 0.9541395)
  expect_equal(res_canned, res_custom)
})
