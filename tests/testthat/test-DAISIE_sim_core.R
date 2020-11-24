context("DAISIE_sim_core")

test_that("new and v1.4a should give same results", {

  tol <- 1e-13
  sim_time <- 10
  n_mainland_species <- 1
  clado_rate <- 1.0
  ext_rate <- 0.5
  carr_cap <- 10
  imm_rate <- 1.0
  ana_rate <- 1.0
  pars <- c(clado_rate, ext_rate, carr_cap, imm_rate, ana_rate)
  rng_seed <- 42
  set.seed(rng_seed)
  new <- DAISIE:::DAISIE_sim_core_constant_rate(
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
  set.seed(rng_seed)
  old <- DAISIE:::DAISIE_sim_core_1_4a(
    time = sim_time,
    mainland_n = n_mainland_species,
    pars = pars
  )

  testthat::expect_true(all(names(new) == names(old)))
  # stt_table has different content
  testthat::expect_true(nrow(new$stt_table) == nrow(old$stt_table))
  # different branching times
  testthat::expect_equal(length(new$branching_times), length(old$branching_times))
  testthat::expect_true(new$stac == old$stac)
  testthat::expect_true(new$missing_species == old$missing_species)
  testthat::expect_true(length(new$all_colonisations) == length(old$all_colonisations))
  testthat::expect_true(new$all_colonisations[[1]]$species_type == old$all_colonisations[[1]]$species_type)

  testthat::expect_true(all(abs(new$stt_table - old$stt_table) < tol))
  testthat::expect_true(all(abs(new$branching_times - old$branching_times) < tol))
  testthat::expect_true(new$all_colonisations[[1]]$brts_miss == old$all_colonisations[[1]]$brts_miss)

  # Frog example
  rng_seed <- 1234
  set.seed(rng_seed)
  time <- 30
  M <- 300
  parsCS <- c(0.437010183, 0.112633464, 36.43883246, 0.00073485, 0)
  new <- DAISIE:::DAISIE_sim_core_constant_rate(
    time = time,
    mainland_n = M,
    pars = parsCS,
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
  set.seed(rng_seed)
  old <- DAISIE:::DAISIE_sim_core_1_4a(
    time = time,
    mainland_n = M,
    pars = parsCS
  )

  testthat::expect_true(all(names(new) == names(old)))
  # stt_table has different content
  testthat::expect_true(nrow(new$stt_table) == nrow(old$stt_table))
  # different branching times
  testthat::expect_equal(length(new$branching_times), length(old$branching_times))
  testthat::expect_true(all(abs(new$stt_table - old$stt_table) < tol))

  for(i in 1:2){
    testthat::expect_true(new$taxon_list[[i]]$stac == old$taxon_list[[i]]$stac)
    testthat::expect_true(new$taxon_list[[i]]$missing_species == old$taxon_list[[i]]$missing_species)
    testthat::expect_true(length(new$taxon_list[[i]]$all_colonisations) == length(old$taxon_list[[i]]$all_colonisations))
    testthat::expect_true(all(abs(new$taxon_list[[i]]$branching_times - old$taxon_list[[i]]$branching_times) < tol))
  }
})

test_that("new and v1.5 should give same results", {

  # Frog example
  tol <- 1e-14
  rng_seed <- 1234
  set.seed(rng_seed)
  time <- 30
  M <- 300
  parsCS <- c(0.437010183, 0.112633464, 36.43883246, 0.00073485, 0)
  new <- DAISIE:::DAISIE_sim_core_constant_rate(
    time = time,
    mainland_n = M,
    pars = parsCS,
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
  set.seed(rng_seed)
  old <- DAISIE:::DAISIE_sim_core_1_5(
    time = time,
    mainland_n = M,
    pars = parsCS
  )
  testthat::expect_true(all(names(new) == names(old)))
  # stt_table has different content
  testthat::expect_true(nrow(new$stt_table) == nrow(old$stt_table))
  # different branching times
  testthat::expect_equal(length(new$branching_times), length(old$branching_times))
  testthat::expect_true(all(abs(new$stt_table - old$stt_table) < tol))

})

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
    DAISIE:::DAISIE_sim_core_constant_rate(
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
  testthat::expect_silent(DAISIE:::DAISIE_sim_core_time_dependent(
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

  ontogeny_sim <- DAISIE_sim_time_dependent(
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
