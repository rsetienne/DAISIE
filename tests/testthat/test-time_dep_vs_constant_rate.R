context("DAISIE_sim_time_dependent first line")

test_that("constant rate output matches time dependent code", {

  # Note: Since both algorithms do not call the RNG an equal number of times,
  # the output of Gillespie runs must necessarily be different. To at least
  # test some of the output, we verify if the first events match, as the first
  # event of each algorithm is sampled from the same number of RNG calls.



# Constant rate code ------------------------------------------------------
  sim_time <- 10
  n_mainland_species <- 1
  clado_rate <- 1.0
  ext_rate <- 0.5
  carr_cap <- 10
  imm_rate <- 1.0
  ana_rate <- 1.0
  pars <- c(clado_rate, ext_rate, carr_cap, imm_rate, ana_rate)
  area_pars <- DAISIE::create_area_pars(
    max_area = 1,
    current_area = 1,
    proportional_peak_t = 0,
    total_island_age = 10,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0
  )
  hyper_pars <- DAISIE::create_hyper_pars(d = 0, x = 0)
  nonoceanic_pars <- c(0, 0)
  rng_seed <- 42
  set.seed(rng_seed)
  constant_rate_out <- DAISIE:::DAISIE_sim_constant_rate(
    time = sim_time,
    M = n_mainland_species,
    pars = pars,
    replicates = 1,
    plot_sims = FALSE,
    verbose = FALSE,
    sample_freq = Inf
  )


  #   Ontogeny code running constant case -----------------------------------
  # We must use the core function to avoid calling calc_peak with a constant
  # area.
  hyper_pars <- DAISIE::create_hyper_pars(d = 0, x = 0)
  set.seed(rng_seed)
  time_dependent_out <- DAISIE:::DAISIE_sim_core_time_dependent(
    time = sim_time,
    mainland_n = n_mainland_species,
    pars = c(clado_rate, ext_rate, carr_cap, imm_rate, ana_rate),
    island_ontogeny = "beta",
    sea_level = "const",
    area_pars = area_pars,
    hyper_pars = hyper_pars,
    nonoceanic_pars = nonoceanic_pars,
    Amax = 1,
    Amin = 1,
    peak = 1
  )

  expect_equal(
    time_dependent_out[[1]][2, ],
    constant_rate_out[[1]][[1]]$stt_all[2, 1:4]
  )

  # Following lines will necessarily be different, see note.
  expect_true(
    !all(time_dependent_out[[1]][3, ] ==
         constant_rate_out[[1]][[1]]$stt_all[3, 1:4])
  )
})

