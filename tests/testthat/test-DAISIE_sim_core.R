context("DAISIE_sim_core")

test_that("new and v1.4 should give same results", {
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
  new <- DAISIE:::DAISIE_sim_core(
    time = sim_time,
    mainland_n = n_mainland_species,
    pars = pars
  )
  set.seed(rng_seed)
  old <- DAISIE:::DAISIE_sim_core_1_4(
    time = sim_time,
    mainland_n = n_mainland_species,
    pars = pars
  )
  #new has init_nonend_spec and init_end_spec in names(new[6:7])
  expect_true(all(names(new[1:5]) == names(old)))
  # stt_table has different content
  expect_true(nrow(new$stt_table) == nrow(old$stt_table))
  # different branching times
  expect_equal(length(new$branching_times), length(old$branching_times))
  expect_true(new$stac == old$stac)
  expect_true(new$missing_species == old$missing_species)
  expect_true(length(new$other_clades_same_ancestor) ==
                length(old$other_clades_same_ancestor))
  expect_true(new$other_clades_same_ancestor[[1]]$species_type ==
                old$other_clades_same_ancestor[[1]]$species_type)

  expect_true(all(new$stt_table == old$stt_table))
  expect_true(all(new$branching_times == old$branching_times))
  expect_true(new$other_clades_same_ancestor[[1]]$brts_miss ==
                old$other_clades_same_ancestor[[1]]$brts_miss)
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
  ddmodel_sim <- 11
  island_type <- "oceanic"
  expect_silent(
    DAISIE:::DAISIE_sim_core(
      time = sim_time,
      mainland_n = n_mainland_species,
      pars = c(clado_rate, ext_rate, carr_cap, imm_rate, ana_rate),
      ddmodel_sim = ddmodel_sim,
      island_type = island_type
    )
  )

})

test_that("Ontogeny oceanic should run silent", {
  set.seed(234567890)
  DAISIE:::DAISIE_sim_core(
    time = 10,
    mainland_n = 1000,
    pars = c(0.0001, 2.2, 0.005, 0.001, 1),
    ddmodel_sim = 11,
    island_type = "oceanic",
    area_pars = create_area_pars(
      max_area = 5000,
      proportional_peak_t = 0.5,
      peak_sharpness = 1,
      total_island_age = 15,
      sea_level_amplitude = 0,
      sea_level_frequency = 0
    ),
    ext_pars = c(1, 100),
    island_ontogeny = "beta",
    sea_level = "const"
  )
  expect_silent(
    DAISIE:::DAISIE_sim_core(
      time = 10,
      mainland_n = 1,
      pars = c(2.5, 2.2, 10, 0.009, 1.01),
      ddmodel_sim = 11,
      island_type = "oceanic",
      area_pars = create_area_pars(5000, 0.2, 1, 15, 0, 0),
      ext_pars = c(1.7, 100),
      island_ontogeny = "beta",
      sea_level = "const"
    )
  )
})

test_that("all species extinct if island dead", {
  ontogeny_sim <- DAISIE:::DAISIE_sim_core(
    time = 10,
    mainland_n = 1000,
    pars = c(0.0001, 2.2, 0.005, 0.001, 1),
    ddmodel_sim = 11,
    island_type = "oceanic",
    area_pars = create_area_pars(
      max_area = 5000,
      proportional_peak_t = 0.5,
      peak_sharpness = 1,
      total_island_age = 10,
      sea_level_amplitude = 0,
      sea_level_frequency = 0
    ),
    ext_pars = c(1, 100),
    island_ontogeny = "beta",
    sea_level = "const"
  )
  last_entry <- ontogeny_sim$stt_table[nrow(ontogeny_sim$stt_table), ]
  expect_true(last_entry[1] == 0)
  expect_true(last_entry[2] == 0)
  expect_true(last_entry[3] == 0)
  expect_true(last_entry[4] == 0)
})

test_that("A non-oceanic run with non-zero sampling should have native
          species on the island", {
            nonoceanic_sim <- DAISIE:::DAISIE_sim_core(
              time = 0.4,
              mainland_n = 1000,
              pars = c(
                2.550687345,
                2.683454548,
                10.0,
                0.00933207,
                1.010073119),
              ddmodel_sim = 11,
              island_type = "nonoceanic",
              nonoceanic_pars = c(0.1, 0.9))
            expect_gt(nonoceanic_sim$stt_table[1, 2], 0)
            expect_gt(nonoceanic_sim$stt_table[1, 3], 0)
          })


test_that("DAISIE_sim_core output is correct", {
  time <- 1
  mainland_n <- 100
  set.seed(17)
  sim_core <- DAISIE_sim_core(time = time,
                              mainland_n = mainland_n,
                              pars = c(2, 2, 20, 0.1, 1))
  expect_true(is.matrix(sim_core$stt_table))
  expect_true(sim_core$stt_table[1, 1] == time)
  expect_true(sim_core$stt_table[nrow(sim_core$stt_table), 1] == 0)
  expect_true(is.numeric(sim_core$taxon_list[[1]]$branching_times))
  expect_true(is.numeric(sim_core$taxon_list[[1]]$stac))
  expect_true(is.numeric(sim_core$taxon_list[[1]]$missing_species))
  expect_true(length(sim_core$taxon_list) == 4)
  expect_true("branching_times" %in% names(sim_core$taxon_list[[1]]))
  expect_true("stac" %in% names(sim_core$taxon_list[[1]]))
  expect_true("missing_species" %in% names(sim_core$taxon_list[[1]]))
  expect_true("init_nonend_spec" %in% names(sim_core$taxon_list[[1]]))
  expect_true("init_end_spec" %in% names(sim_core$taxon_list[[1]]))
  expect_true("carrying_capacity" %in% names(sim_core$taxon_list[[1]]))
})

test_that("DAISIE_sim_core fails when pars[4] == 0 &&
          island_type == 'oceanic'", {
            expect_error(DAISIE_sim_core(time = 1,
                                         mainland_n = 100,
                                         pars = c(2, 2, 20, 0, 1),
                                         island_type = "oceanic"))
          })

test_that("!is.null(area_pars) && island_ontogeny == 'const'", {
  expect_error(DAISIE_sim_core(time = 1,
                               mainland_n = 100,
                               pars = c(2, 2, 20, 0, 1),
                               island_ontogeny = "const",
                               area_pars = create_area_pars(max_area = 1,
                                                        proportional_peak_t = 1,
                                                        peak_sharpness = 1,
                                                        total_island_age = 1)))
})

test_that("split-rate model runs silent and
          gives correct output", {
  expect_silent(DAISIE_sim_core(time = 1,
                                mainland_n = 1,
                                pars = c(1,1,1,1,1,1,1,1,1,1),
                                pars_shift = TRUE,
                                shift_times = 5))
})
