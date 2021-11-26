test_that("Clean run should be silent", {

  set.seed(42)
  n_mainland_species <- 1
  sim_time <- 10
  clado_rate <- 1.0
  ext_rate <- 0.1
  carr_cap <- 4
  imm_rate <- 1.0
  ana_rate <- 1.0

  testthat::expect_silent(
    DAISIE:::DAISIE_sim_core_cr(
      time = sim_time,
      mainland_n = n_mainland_species,
      pars = c(clado_rate, ext_rate, carr_cap, imm_rate, ana_rate),
      nonoceanic_pars = c(0, 0),
      hyper_pars = create_hyper_pars(d = 0, x = 0),
      area_pars = create_area_pars(
        max_area = 1,
        current_area = 1,
        proportional_peak_t = 0,
        total_island_age = 0,
        sea_level_amplitude = 0,
        sea_level_frequency = 0,
        island_gradient_angle = 0
      )
    )
  )

})

test_that("Ontogeny simulation should run silent", {
  set.seed(234567890)
  testthat::expect_silent(DAISIE:::DAISIE_sim_core_time_dep(
    time = 10,
    mainland_n = 1000,
    pars = c(0.0001, 2.2, 0.005, 0.001, 1),
    nonoceanic_pars = c(0, 0),
    area_pars = create_area_pars(
      max_area = 5000,
      current_area = 2500,
      proportional_peak_t = 0.5,
      total_island_age = 15,
      sea_level_amplitude = 0,
      sea_level_frequency = 0,
      island_gradient_angle = 0),
    island_ontogeny = "beta",
    Amax = 5000,
    Amin = 0,
    peak = 1,
    hyper_pars = create_hyper_pars(
      d = 0.2,
      x = 0.15
    )
  )
  )
})

test_that("all species extinct if island dead", {
  set.seed(234567890)

  ontogeny_sim <- DAISIE_sim_time_dep(
    time = 10,
    M = 1,
    pars = c(0.0001, 2.5, 0.005, 0.01, 1),
    nonoceanic_pars = c(0, 0),
    replicates = 1,
    area_pars = create_area_pars(
      max_area = 5000,
      current_area = 0.0000001,
      proportional_peak_t = 0.5,
      total_island_age = 10.000000001,
      sea_level_amplitude = 0,
      sea_level_frequency = 0,
      island_gradient_angle = 0),
    island_ontogeny = "beta",
    hyper_pars = create_hyper_pars(
      d = 0.2,
      x = 0.15
    ),
    verbose = FALSE
  )
  last_entry <- ontogeny_sim[[1]][[1]]$stt_all[nrow(ontogeny_sim[[1]][[1]]$stt_all), ]
  testthat::expect_true(last_entry[1] == 0)
  testthat::expect_true(last_entry[2] == 0)
  testthat::expect_true(last_entry[3] == 0)
  testthat::expect_true(last_entry[4] == 0)
})
