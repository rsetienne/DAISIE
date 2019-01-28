context("DAISIE_calc_clade_imm_rate")

test_that("use, below carrying capacity", {

  ps_imm_rate <- 12.34
  n_island_species <- 42
  n_mainland_species <- 21
  carr_cap <- 314
  created <- DAISIE_calc_clade_imm_rate(
    ps_imm_rate = ps_imm_rate, 
    n_island_species = n_island_species, 
    n_mainland_species = n_mainland_species, 
    carr_cap = carr_cap
  )
  expected <- n_mainland_species * ps_imm_rate * (1.0 - (n_island_species / carr_cap))
  expect_equal(created, expected)
})

test_that("use, above carrying capacity", {

  ps_imm_rate <- 12.34
  n_mainland_species <- 21
  carr_cap <- 314
  n_island_species <- carr_cap + 42
  n_mainland_species <- carr_cap + 21
  created <- DAISIE_calc_clade_imm_rate(
    ps_imm_rate = ps_imm_rate, 
    n_island_species = n_island_species, 
    n_mainland_species = n_mainland_species, 
    carr_cap = carr_cap
  )
  expected <- 0.0
  expect_equal(created, expected)
})
