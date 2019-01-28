context("DAISIE_sim_with_mainland")

test_that("Without mainland processes, should do the same as classic model", {
  # Create parameters
  sim_time <- 10
  n_mainland_species <- 10
  clado_rate <- 1.0
  ext_rate <- 0.1
  carr_cap <- 10
  imm_rate <- 0.1
  ana_rate <- 1.0
  n_replicates <- 3
  seed <- 42
  testit::assert(sim_time > 0.0)
  testit::assert(n_mainland_species > 0)
  testit::assert(clado_rate >= 0.0)
  testit::assert(ext_rate >= 0.0)
  testit::assert(carr_cap > 0)
  testit::assert(imm_rate > 0.0)
  testit::assert(ana_rate >= 0.0)
  testit::assert(n_replicates > 0)
  pars <- c(clado_rate, ext_rate, carr_cap, imm_rate, ana_rate)
  
  # Run
  set.seed(seed)
  expected <- DAISIE_sim(
    time = sim_time,
    M = n_mainland_species,
    pars = pars,
    replicates = n_replicates,
    verbose = FALSE,
    plot_sims = FALSE
  )
  set.seed(seed)
  created <- DAISIE_sim_with_mainland(
    time = sim_time,
    M = n_mainland_species,
    pars = pars,
    replicates = n_replicates
  )
  
  # Compare for being exactly the same
  # This test will soften up in the future.
  # In the future, the STT's are compared to be equivalent
  testit::assert(n_replicates == length(created))
  expect_equal(length(expected), length(created))
  for (replicate_index in seq_along(expected)) {
    this_expected <- expected[[replicate_index]]
    this_created <- created[[replicate_index]]
    expect_equal(length(this_expected), length(this_created))
    for (element_index in seq_along(this_expected)) 
    {
      expected_element <- this_expected[[element_index]]
      created_element <- this_created[[element_index]]
      expect_equal(expected_element, created_element)
    }
  }
})

test_that("can plot", {
  # Create parameters
  sim_time <- 10
  n_mainland_species <- 10
  clado_rate <- 1.0
  ext_rate <- 0.1
  carr_cap <- 10
  imm_rate <- 1.0
  ana_rate <- 1.0
  n_replicates <- 3
  seed <- 42
  testit::assert(sim_time > 0.0)
  testit::assert(n_mainland_species > 0)
  testit::assert(clado_rate >= 0.0)
  testit::assert(ext_rate >= 0.0)
  testit::assert(carr_cap > 0)
  testit::assert(imm_rate > 0.0)
  testit::assert(ana_rate >= 0.0)
  testit::assert(n_replicates > 0)
  pars <- c(clado_rate, ext_rate, carr_cap, imm_rate, ana_rate)
  
  set.seed(seed)
  results <- DAISIE_sim_with_mainland(
    time = sim_time,
    M = n_mainland_species,
    pars = pars,
    replicates = n_replicates
  )
  DAISIE_plot_sims(
    island_replicates = results, 
    use_dev_new = FALSE, 
    plot_plus_one = FALSE
  )
})
