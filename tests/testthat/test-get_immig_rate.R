context("get_immig_rate")

test_that("classic behavior", {
  
  carr_cap <- 10
  ps_imm_rate <- 0.1
  n_island_species <- 5
  n_mainland_species <- 2
  expected <- DAISIE_calc_clade_imm_rate(
    ps_imm_rate = ps_imm_rate,
    n_island_species = n_island_species,
    n_mainland_species = n_mainland_species,
    carr_cap = carr_cap
  )
  created <- get_immig_rate(
    timeval = 1.0,
    totaltime = 10.0,
    gam = ps_imm_rate,
    Apars =  NULL,
    island_ontogeny = 0,
    island_spec = matrix(data = NA, nrow = n_island_species, ncol = 1),
    K = carr_cap,
    mainland_n = n_mainland_species
  )
  expect_equal(expected, created)
})
