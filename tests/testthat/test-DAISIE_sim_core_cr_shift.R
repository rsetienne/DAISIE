test_that("split-rate model runs silent and
          gives correct output", {
            set.seed(1)
            area_pars <- create_area_pars(
              max_area = 1,
              current_area = 1,
              proportional_peak_t = 0,
              total_island_age = 0,
              sea_level_amplitude = 0,
              sea_level_frequency = 0,
              island_gradient_angle = 0)
            nonoceanic_pars <- c(0, 0)
            hyper_pars <- create_hyper_pars(d = 0, x = 0)
            testthat::expect_silent(
              DAISIE_sim_core_cr_shift(
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
  area_pars <- create_area_pars(
    max_area = 1,
    current_area = 1,
    proportional_peak_t = 0,
    total_island_age = 0,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0)
  nonoceanic_pars <- c(0, 0)
  hyper_pars <- create_hyper_pars(d = 0, x = 0)
  testthat::expect_error(DAISIE_sim_core_cr_shift(
    time = 1,
    mainland_n = 1,
    pars = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    shift_times = 5,
    area_pars = area_pars)
  )
})

test_that("abuse split-rate model with gamma = 0", {
  area_pars <- create_area_pars(
    max_area = 1,
    current_area = 1,
    proportional_peak_t = 0,
    total_island_age = 0,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0)
  nonoceanic_pars <- c(0, 0)
  hyper_pars <- create_hyper_pars(d = 0, x = 0)
  testthat::expect_error(DAISIE_sim_core_cr_shift(
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


test_that("a shift at the island origin starts and stays in the second rate set", {
  area_pars <- create_area_pars(
    max_area = 1,
    current_area = 1,
    proportional_peak_t = 0,
    total_island_age = 0,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0)
  hyper_pars <- create_hyper_pars(d = 0, x = 0)
  total_time <- 1
  mainland_n <- 5
  # Rate set with (practically) no immigration, and one with immigration so high
  # that the island is colonised almost immediately
  no_immigration <- c(1e-9, 1e-9, 20, 1e-9, 1e-9)
  high_immigration <- c(1e-9, 1e-9, 20, 20, 1e-9)
  n_species_at_present <- function(island) {
    stt_table <- island$stt_table
    sum(stt_table[nrow(stt_table), c("nI", "nA", "nC")])
  }

  # shift_times = total_time is a shift at time 0: the island must be in the
  # second rate set for the whole simulation
  set.seed(1)
  island <- DAISIE_sim_core_cr_shift(
    time = total_time,
    mainland_n = mainland_n,
    pars = c(no_immigration, high_immigration),
    shift_times = total_time,
    area_pars = area_pars,
    hyper_pars = hyper_pars
  )
  testthat::expect_gt(n_species_at_present(island), 0)

  set.seed(1)
  island <- DAISIE_sim_core_cr_shift(
    time = total_time,
    mainland_n = mainland_n,
    pars = c(high_immigration, no_immigration),
    shift_times = total_time,
    area_pars = area_pars,
    hyper_pars = hyper_pars
  )
  testthat::expect_equal(n_species_at_present(island), 0)
})
