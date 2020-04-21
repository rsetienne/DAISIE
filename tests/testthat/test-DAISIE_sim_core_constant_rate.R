context("DAISIE_sim_core_constant_rate")

test_that("Clean run should be silent", {
  set.seed(42)
  n_mainland_species <- 1
  sim_time <- 10
  clado_rate <- 1.0
  ext_rate <- 0.1
  carr_cap <- 4
  imm_rate <- 1.0
  ana_rate <- 1.0
  area_pars <- DAISIE::create_area_pars(
    max_area = 1,
    current_area = 1,
    proportional_peak_t = 0,
    total_island_age = 0,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0)
  hyper_pars <- create_hyper_pars(d = 0, x = 0)
  nonoceanic_pars <- c(0, 0)
  expect_silent(
    DAISIE:::DAISIE_sim_core_constant_rate(
      time = sim_time,
      mainland_n = n_mainland_species,
      pars = c(clado_rate, ext_rate, carr_cap, imm_rate, ana_rate),
      area_pars = area_pars,
      hyper_pars = hyper_pars,
      nonoceanic_pars = nonoceanic_pars
    )
  )
})

test_that("A non-oceanic run with non-zero sampling should have native
          species on the island", {
            area_pars <- DAISIE::create_area_pars(
              max_area = 1,
              current_area = 1,
              proportional_peak_t = 0,
              total_island_age = 0,
              sea_level_amplitude = 0,
              sea_level_frequency = 0,
              island_gradient_angle = 0)
            hyper_pars <- create_hyper_pars(d = 0, x = 0)
            nonoceanic_sim <- DAISIE:::DAISIE_sim_core_constant_rate(
              time = 0.4,
              mainland_n = 1000,
              pars = c(
                2.550687345,
                2.683454548,
                10.0,
                0.00933207,
                1.010073119),
              area_pars = area_pars,
              hyper_pars = hyper_pars,
              nonoceanic_pars = c(0.1, 0.9)
            )
            expect_gt(nonoceanic_sim$stt_table[1, 2], 0)
            expect_gt(nonoceanic_sim$stt_table[1, 3], 0)
})


test_that("DAISIE_sim_core output is correct", {
  time <- 1
  mainland_n <- 100
  set.seed(5)
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
  sim_core <- DAISIE:::DAISIE_sim_core_constant_rate(
    time = time,
    mainland_n = mainland_n,
    pars = c(2, 2, 20, 0.1, 1),
    area_pars = area_pars,
    hyper_pars = hyper_pars,
    nonoceanic_pars = nonoceanic_pars
    )
  expect_true(is.matrix(sim_core$stt_table))
  expect_true(sim_core$stt_table[1, 1] == time)
  expect_true(sim_core$stt_table[nrow(sim_core$stt_table), 1] == 0)
  expect_true(is.numeric(sim_core$taxon_list[[1]]$branching_times))
  expect_true(is.numeric(sim_core$taxon_list[[1]]$stac))
  expect_true(is.numeric(sim_core$taxon_list[[1]]$missing_species))
  expect_true(length(sim_core$taxon_list) == 2)
  expect_true("branching_times" %in% names(sim_core$taxon_list[[1]]))
  expect_true("stac" %in% names(sim_core$taxon_list[[1]]))
  expect_true("missing_species" %in% names(sim_core$taxon_list[[1]]))
})

test_that("DAISIE_sim_core with land-bridge starting at time = 0 for CS uses
          the second parameter set at time = 0", {
            area_pars <- DAISIE::create_area_pars(
              max_area = 1,
              current_area = 1,
              proportional_peak_t = 0,
              total_island_age = 0,
              sea_level_amplitude = 0,
              sea_level_frequency = 0,
              island_gradient_angle = 0)
            hyper_pars <- create_hyper_pars(d = 0, x = 0)
            expect_silent(
              DAISIE:::DAISIE_sim_core_constant_rate_shift(
                time = 10,
                mainland_n = 1,
                pars = c(1, 1, 10, 0.1, 1, 2, 2, 20, 0.2, 1),
                shift_times = 10,
                area_pars = area_pars,
                hyper_pars = hyper_pars
              )
            )
          })

test_that("DAISIE_sim_core fails when pars[4] == 0 &&
          nonoceanic_pars[1] == 0", {
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
            expect_error(
              DAISIE:::DAISIE_sim_core_constant_rate(
                time = 1,
                mainland_n = 100,
                pars = c(2, 2, 20, 0, 1),
                area_pars = area_pars,
                hyper_pars = hyper_pars,
                nonoceanic_pars = nonoceanic_pars
              )
            )
})

