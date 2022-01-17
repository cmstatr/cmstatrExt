test_that("equiv_sample", {
  k <- k_equiv_two_sample(18, 5, 0.05)
  expect_equal(k[1], 2.867903, tolerance = 1e-3)
  expect_equal(k[2], 1.019985, tolerance = 1e-3)

  # the following used to fail
  k_equiv_two_sample(100, 6, 0.05)
})
