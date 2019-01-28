context("DAISIE_sim_core_checked")

test_that("use", {

  set.seed(42)
  n_mainland_species <- 1
  sim_time <- 10
  result <- NULL
  expect_silent(
    result <- DAISIE::DAISIE_sim_core_checked(
      sim_time = sim_time, 
      n_mainland_species = n_mainland_species, 
      clado_rate = 1.0, 
      ext_rate = 0.1,
      carr_cap = 4,
      imm_rate = 1.0,
      ana_rate = 1.0
    )
  )
  expect_true("stt_table" %in% names(result))
  expect_true("branching_times" %in% names(result))
  expect_true("stac" %in% names(result))
  expect_true("missing_species" %in% names(result))
  if (1 == 2) {
    island_replicates <- list()
    island_replicates[[1]] <- list()
    island_replicates[[1]][[1]] <- result
    # May also return 'other_clades_same_ancestor'
    island_replicates <- DAISIE_format_CS(
      island_replicates = island_replicates,
      time = sim_time,
      M = n_mainland_species,
      sample_freq = 25,
      verbose = FALSE
    )
    DAISIE_plot_sims(island_replicates, use_dev_new = FALSE)
    DAISIE_plot_island(island = island_replicates[[1]], island_age = sim_time) 
  }
})
