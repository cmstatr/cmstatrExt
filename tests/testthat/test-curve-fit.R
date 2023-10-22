suppressMessages(library(dplyr))

test_that("average_curve_lm produces expected output", {
  res <- pa12_tension %>%
    average_curve_lm(
      Coupon,
      Stress ~ I(Strain) + I(Strain^2) + I(Strain^3) + 0,
      n_bins = 100
    )

  expect_equal(res$n_bins, 100)
  expect_equal(as.character(res$y_var), "Stress")
  expect_equal(as.character(res$x_var), "Strain")
  expect_equal(res$data, pa12_tension)
  expect_equal(nrow(res$binned_data), 400)  # 4 groups
  expect_snapshot(print(res))
  expect_snapshot(summary(res))

  augmented_dat <- augment(res)  # newdata is NULL, extrapolate is FALSE
  expect_equal(nrow(augmented_dat), nrow(pa12_tension))
  expect_length(augmented_dat, 6)  # original 3 columns and 3 more
  expect_lte(
    (augmented_dat %>%
       filter(!is.na(`.fit`)) %>%
       summarise(max_strain = max(Strain)) %>%
       select(c(`max_strain`)))[[1]],
    res$max_x
  )
  expect_equal(
    augmented_dat %>%
      filter(Strain <= res$max_x) %>%
      nrow(),
    augmented_dat %>%
      filter(!is.na(`.fit`)) %>%
      nrow()
  )

  augmented_dat <- augment(res, extrapolate = TRUE)  # newdata is NULL
  expect_equal(nrow(augmented_dat), nrow(pa12_tension))
  expect_length(augmented_dat, 6)  # original 3 columns and 3 more
  expect_equal(
    (augmented_dat %>%
       filter(!is.na(`.fit`)) %>%
       summarise(max_strain = max(Strain)) %>%
       select(c(`max_strain`)))[[1]],
    (pa12_tension %>%
       summarise(max_strain = max(Strain)) %>%
       select(c(`max_strain`)))[[1]]
  )
  expect_equal(
    augmented_dat %>%
      filter(!is.na(`.fit`)) %>%
      nrow(),
    pa12_tension %>%
      nrow()
  )

  dat <- pa12_tension %>%
    filter(`Coupon` == "Coupon 1")
  augmented_dat <- augment(res, newdata = dat, extrapolate = TRUE)
  expect_equal(nrow(augmented_dat), nrow(dat))
  expect_length(augmented_dat, 6)  # original 3 columns and 3 more
  expect_equal(
    (augmented_dat %>%
       filter(!is.na(`.fit`)) %>%
       summarise(max_strain = max(Strain)) %>%
       select(c(`max_strain`)))[[1]],
    (dat %>%
       summarise(max_strain = max(Strain)) %>%
       select(c(`max_strain`)))[[1]]
  )
  expect_equal(
    augmented_dat %>%
      filter(!is.na(`.fit`)) %>%
      nrow(),
    dat %>%
      nrow()
  )
})

test_that("average_curve_lm produces expected errors and warnings", {
  expect_error(
    average_curve_lm(pa12_tension, Sample, ~ Strain),
    "LHS.+one variable"
  )
  expect_error(
    average_curve_lm(pa12_tension, Sample, Stress + Strain ~ Strain),
    "LHS.+one variable"
  )

  expect_error(
    average_curve_lm(pa12_tension, Sample, Strain ~ NULL),
    "RHS.+one variable"
  )
  expect_error(
    average_curve_lm(pa12_tension, Sample, Strain ~ Stress + Strain),
    "RHS.+one variable"
  )

  expect_error(
    average_curve_lm(pa12_tension, Coupon, Strain ~ Strain),
    "different variable"
  )

  expect_warning(
    pa12_tension %>%
      filter(Strain < 0.1 | Strain > 0.12) %>%
      average_curve_lm(Coupon, Stress ~ Strain),
    "empty"
  )

  expect_warning(
    pa12_tension %>%
      filter(!(Coupon == "Coupon 1" & (Strain > 0.1 & Strain < 0.12))) %>%
      average_curve_lm(Coupon, Stress ~ Strain),
    "empty"
  )

  expect_error(
    pa12_tension %>%
      mutate(Strain = -Strain) %>%
      average_curve_lm(Coupon, Stress ~ Strain),
    "No positive"
  )

  expect_warning(
    pa12_tension %>%
      mutate(Strain = Strain - 0.001) %>%
      average_curve_lm(Coupon, Stress ~ Strain),
    "ignored"
  )
})

test_that("average_curve_optim produces expected output", {
  res_opt <- average_curve_optim(
    pa12_tension,
    Coupon,
    Strain,
    Stress,
    function(strain, par) {
      sum(par * c(strain, strain^2, strain^3))
    },
    c(c1 = 1, c2 = 1, c3 = 1),  # TODO: Document this
    n_bins = 100
  )

  expect_equal(res_opt$n_bins, 100)
  expect_equal(as.character(res_opt$y_var), "Stress")
  expect_equal(as.character(res_opt$x_var), "Strain")
  expect_equal(res_opt$data, pa12_tension)
  expect_equal(nrow(res_opt$binned_data), 400)  # 4 groups
  expect_snapshot(print(res_opt))
})

test_that("average_curve_optim and average_curve_lm produce similar results", {
  res_opt <- average_curve_optim(
    pa12_tension,
    Coupon,
    Strain,
    Stress,
    function(strain, par) {
      sum(par * c(strain, strain^2, strain^3))
    },
    c(c1 = 1, c2 = 1, c3 = 1),
    n_bins = 100
  )

  res_lm <- average_curve_lm(
    pa12_tension,
    Coupon,
    Stress ~ I(Strain) + I(Strain^2) + I(Strain^3) + 0,
    n_bins = 100
  )

  expect_equal(res_opt$fit_optim$par, res_lm$fit_lm$coefficients,
               ignore_attr = TRUE, tolerance = 0.0001)

  expect_equal(augment(res_opt), augment(res_lm),
               ignore_attr = TRUE, tolerance = 0.0001)
})
