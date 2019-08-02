context("DAISIE_half_life")
skip("WIP")
test_that("DAISIE_half_life calculates half-life", {
  sim_time <- 10
  n_mainland_species <- 1000
  clado_rate <- 1.0
  ext_rate <- 0.5
  carr_cap <- 10
  imm_rate <- 1.0
  ana_rate <- 1.0
  pars <- c(clado_rate, ext_rate, carr_cap, imm_rate, ana_rate)
  sim_core_half_life <- DAISIE_sim_core(time = 3,
                                        mainland_n = 1000,
                                        pars = pars)
  sim_core_half_life <- DAISIE_half_life(sim_core = sim_core_half_life)
  expect_true(is.numeric(sim_core_half_life))
  expect_true(sim_core_half_life != 0)
})
skip("WIP")
test_that("DAISIE_half_life changes for diversity-dependent extinction", {
  sim_time <- 10
  n_mainland_species <- 1000
  clado_rate <- 1.0
  ext_rate <- 0.5
  carr_cap <- 10
  imm_rate <- 1.0
  ana_rate <- 1.0
  pars <- c(clado_rate, ext_rate, carr_cap, imm_rate, ana_rate)
  divdep_DImu <- c("lac", "gam")
  divdep_DDmu <- c("lac", "mu", "gam")
  sim_core_DImu <- DAISIE_sim_core(time = 3,
                                    mainland_n = 1000,
                                    pars = pars,
                                    divdep = divdep_DImu)
  half_life_DImu <- DAISIE_half_life(sim_core = sim_core_DImu)
  sim_core_DDmu <- DAISIE_sim_core(time = 3,
                                    mainland_n = 1000,
                                    pars = pars,
                                    divdep = divdep_DDmu)
  
  half_life_DDmu <- DAISIE_half_life(sim_core = sim_core_DDmu)
  expect_true(half_life_DImu != half_life_DDmu)
})
skip("WIP")
test_that("DAISIE_half_life and DAISIE_avg_half_life converge", {
  sim_time <- 10
  n_mainland_species <- 1000
  clado_rate <- 1.0
  ext_rate <- 0.5
  carr_cap <- 10 
  imm_rate <- 1.0
  ana_rate <- 1.0
  pars <- c(clado_rate, ext_rate, carr_cap, imm_rate, ana_rate)
  sim_core_half_life <- DAISIE_sim_core(time = 3,
                                        mainland_n = 1000,
                                        pars = pars)
  sim_core_half_life <- DAISIE_half_life(sim_core = sim_core_half_life)
  expect_true(is.numeric(sim_core_half_life))
  expect_true(sim_core_half_life != 0)
})