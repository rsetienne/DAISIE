context("DAISIE_sim_core_constant_rate_shift")

test_that("split-rate model runs silent and
          gives correct output", {
            set.seed(1)
            area_pars <- DAISIE::create_area_pars(
              max_area = 1,
              current_area = 1,
              proportional_peak_t = 0,
              total_island_age = 0,
              sea_level_amplitude = 0,
              sea_level_frequency = 0,
              island_gradient_angle = 0)
            nonoceanic_pars <- c(0, 0)
            hyper_pars <- create_hyper_pars(d = 0, x = 0)
            expect_silent(
              DAISIE:::DAISIE_sim_core_constant_rate_shift(
                time = 10,
                mainland_n = 1,
                pars = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                shift_times = 5,
                area_pars = area_pars,
                hyper_pars = hyper_pars,
                nonoceanic_pars = nonoceanic_pars
              )
            )
})

test_that("abuse split-rate model with time smaller than shift_times", {
  area_pars <- DAISIE::create_area_pars(
    max_area = 1,
    current_area = 1,
    proportional_peak_t = 0,
    total_island_age = 0,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0)
  nonoceanic_pars <- c(0, 0)
  hyper_pars <- create_hyper_pars(d = 0, x = 0)
  expect_error(DAISIE:::DAISIE_sim_core_constant_rate_shift(
    time = 1,
    mainland_n = 1,
    pars = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    shift_times = 5,
    area_pars = area_pars)
  )
})

test_that("abuse split-rate model with gamma = 0", {
  area_pars <- DAISIE::create_area_pars(
    max_area = 1,
    current_area = 1,
    proportional_peak_t = 0,
    total_island_age = 0,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0)
  nonoceanic_pars <- c(0, 0)
  hyper_pars <- create_hyper_pars(d = 0, x = 0)
  expect_error(DAISIE:::DAISIE_sim_core_constant_rate_shift(
    time = 1,
    mainland_n = 1,
    pars = c(1, 1, 1, 0, 1, 1, 1, 1, 1, 1),
    shift_times = 5,
    area_pars = area_pars),
    regexp =
      "Island has no species and the rate of
    colonisation is zero. Island cannot be colonised."
  )
})

