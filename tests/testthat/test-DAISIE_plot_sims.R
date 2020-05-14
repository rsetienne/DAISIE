context("DAISIE_plot_sims")

test_that("Example 1", {

  data(islands_1type_1000reps)
  DAISIE_plot_sims(
    island_replicates = islands_1type_1000reps
  )
})

test_that("Example 2", {
  data(islands_2types_1000reps)
  DAISIE_plot_sims(
    island_replicates = islands_2types_1000reps
  )
})

test_that("use", {

  set.seed(42)
  n_mainland_species <- 1
  sim_time <- 10
  result <- DAISIE:::DAISIE_sim_core_checked(
    sim_time = sim_time,
    n_mainland_species = n_mainland_species,
    clado_rate = 1.0,
    ext_rate = 0.1,
    carr_cap = 4,
    imm_rate = 1.0,
    ana_rate = 1.0
  )

  island_replicates <- list()
  island_replicates[[1]] <- list()
  island_replicates[[1]][[1]] <- result
  # May also return 'other_clades_same_ancestor'

  island_replicates <- DAISIE:::DAISIE_format_CS(
    island_replicates = island_replicates,
    time = sim_time,
    M = n_mainland_species,
    sample_freq = 25,
    verbose = FALSE
  )
  DAISIE:::DAISIE_plot_sims(
    island_replicates = island_replicates,
    plot_plus_one = FALSE
  )
})

test_that("Plot plus one", {

  data(islands_1type_1000reps)
  DAISIE_plot_sims(
    island_replicates = islands_1type_1000reps,
    plot_plus_one = TRUE
  )
  DAISIE_plot_sims(
    island_replicates = islands_1type_1000reps,
    plot_plus_one = FALSE
  )
})

