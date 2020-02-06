context("DAISIE_sim_core")

test_that("Clean run should be silent", {
  set.seed(42)
  n_mainland_species <- 1
  sim_time <- 10
  clado_rate <- 1.0
  ext_rate <- 0.1
  carr_cap <- 4
  imm_rate <- 1.0
  ana_rate <- 1.0
  expect_silent(
    DAISIE:::DAISIE_sim_core_constant_rate(
      time = sim_time,
      mainland_n = n_mainland_species,
      pars = c(clado_rate, ext_rate, carr_cap, imm_rate, ana_rate)
    )
  )
})

test_that("Ontogeny oceanic should run silent IW", {
  set.seed(234567890)
  expect_silent(
    DAISIE:::DAISIE_sim_core_time_dependent(
      time = 10,
      mainland_n = 100,
      pars = c(0.0001, 2.2, 0.005, 0.001, 1),
      area_pars = create_area_pars(
        max_area = 5000,
        proportional_peak_t = 0.5,
        peak_sharpness = 1,
        total_island_age = 15,
        sea_level_amplitude = 0,
        sea_level_frequency = 0,
        island_gradient_angle = 0
      ),
      ext_pars = c(1, 100),
      island_ontogeny = "beta",
      sea_level = "const",
      extcutoff = 20 #1000
    )
  )
})

test_that("Ontogeny oceanic should run silent CS", {
  set.seed(420)
  expect_silent(
    DAISIE:::DAISIE_sim_core_time_dependent(
      time = 10,
      mainland_n = 1,
      pars = c(0.0001, 2.2, 0.005, 0.001, 1),
      area_pars = create_area_pars(
        max_area = 5000,
        proportional_peak_t = 0.5,
        peak_sharpness = 1,
        total_island_age = 15,
        sea_level_amplitude = 0,
        sea_level_frequency = 0,
        island_gradient_angle = 0
      ),
      ext_pars = c(1, 100),
      island_ontogeny = "beta",
      sea_level = "const"
    )
  )
})

test_that("all species extinct if island dead", {
  ontogeny_sim <- DAISIE:::DAISIE_sim_core_time_dependent(
    time = 10,
    mainland_n = 1000,
    pars = c(0.0001, 2.2, 0.005, 0.001, 1),
    area_pars = create_area_pars(
      max_area = 5000,
      proportional_peak_t = 0.5,
      peak_sharpness = 1,
      total_island_age = 10,
      sea_level_amplitude = 0,
      sea_level_frequency = 0,
      island_gradient_angle = 0
    ),
    ext_pars = c(1, 100),
    island_ontogeny = "beta",
    sea_level = "const",
    extcutoff = 20 #1000
  )
  last_entry <- ontogeny_sim$stt_table[nrow(ontogeny_sim$stt_table), ]
  expect_true(last_entry[1] == 0)
  expect_true(last_entry[2] == 0)
  expect_true(last_entry[3] == 0)
  expect_true(last_entry[4] == 0)
})

test_that("A non-oceanic run with non-zero sampling should have native
          species on the island", {
            nonoceanic_sim <- DAISIE:::DAISIE_sim_core_constant_rate(
              time = 0.4,
              mainland_n = 1000,
              pars = c(
                2.550687345,
                2.683454548,
                10.0,
                0.00933207,
                1.010073119),
              nonoceanic_pars = c(0.1, 0.9))
            expect_gt(nonoceanic_sim$stt_table[1, 2], 0)
            expect_gt(nonoceanic_sim$stt_table[1, 3], 0)
          })


test_that("DAISIE_sim_core output is correct", {
  time <- 1
  mainland_n <- 100
  set.seed(5)
  sim_core <- DAISIE:::DAISIE_sim_core_constant_rate(time = time,
                              mainland_n = mainland_n,
                              pars = c(2, 2, 20, 0.1, 1))
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
  expect_silent(DAISIE:::DAISIE_sim_core_constant_rate_shift(time = 10,
                                mainland_n = 1,
                                pars = c(1, 1, 10, 0.1, 1, 2, 2, 20, 0.2, 1),
                                shift_times = 10))
})

test_that("DAISIE_sim_core fails when pars[4] == 0 &&
          nonoceanic_pars[1] == 0", {
            expect_error(DAISIE:::DAISIE_sim_core_constant_rate(time = 1,
                                         mainland_n = 100,
                                         pars = c(2, 2, 20, 0, 1)))
          })

test_that("!is.null(area_pars) && island_ontogeny == 'const'", {
  expect_error(DAISIE:::DAISIE_sim_core_time_dependent(time = 1,
                               mainland_n = 100,
                               pars = c(2, 2, 20, 0, 1),
                               island_ontogeny = 0,
                               sea_level = 0,
                               area_pars = create_area_pars(
                                 max_area = 1,
                                 proportional_peak_t = 1,
                                 peak_sharpness = 1,
                                 total_island_age = 1,
                                 sea_level_amplitude = 0,
                                 sea_level_frequency = 0,
                                 island_gradient_angle = 0
                               )
  ), regexp = "area_pars specified for constant island_ontogeny and sea_level.
         Run DAISIE_sim_constant_rate instead.")
})

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

test_that("(is.null(ext_pars) || is.null(area_pars)) &&
          (island_ontogeny != 0 || sea_level != 0)", {

            time <- 10
            mainland_n <- 1000
            pars <- c(0.0001, 2.2, 0.005, 0.001, 1)
            area_pars <- NULL
            ext_pars <- NULL
            island_ontogeny <- 1
            sea_level <- 1

            expect_error(
              DAISIE:::DAISIE_sim_core_time_dependent(
                time = time,
                mainland_n = mainland_n,
                pars = pars,
                area_pars = area_pars,
                ext_pars = ext_pars,
                island_ontogeny = island_ontogeny,
                sea_level = sea_level
              ), regexp =
                "Island ontogeny and/or sea level specified but area parameters
    and/or extinction parameters not available. Please either set
    island_ontogeny and sea_level to NULL, or specify area_pars and ext_pars."
            )
})

