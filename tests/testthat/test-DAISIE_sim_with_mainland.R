context("DAISIE_sim_with_mainland")

test_that("Mainland sim call throws error", {
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
  
  expect_error(
    DAISIE:::DAISIE_sim_with_mainland(
      time = sim_time,
      M = n_mainland_species,
      pars = pars,
      replicates = n_replicates,
      mainland_params = pars),
    "Mainland speciation not implemented yet")
})
