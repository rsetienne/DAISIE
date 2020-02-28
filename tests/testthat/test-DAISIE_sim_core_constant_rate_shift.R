context("DAISIE_sim_core_constant_rate_shift")

test_that("split-rate model runs silent and
          gives correct output", {
            set.seed(1)
            expect_silent(DAISIE:::DAISIE_sim_core_constant_rate_shift(time = 10,
                                                                       mainland_n = 1,
                                                                       pars = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                                                                       shift_times = 5))
          })

test_that("abuse split-rate model with time smaller than shift_times", {
  expect_error(DAISIE:::DAISIE_sim_core_constant_rate_shift(time = 1,
                                                            mainland_n = 1,
                                                            pars = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                                                            shift_times = 5))
})

