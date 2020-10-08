context("get_immig_rate")

test_that("use area constant diversity-independent", {
  carr_cap <- Inf
  ps_imm_rate <- 0.1
  n_island_species <- 5
  n_mainland_species <- 1
  hyper_pars <- DAISIE:::create_hyper_pars(0, 0)
  area <- 1
  created <- DAISIE:::get_immig_rate(
    gam = ps_imm_rate,
    A = area,
    num_spec = n_island_species,
    K = carr_cap,
    mainland_n = n_mainland_species)
  expected <- DAISIE_calc_clade_imm_rate(
    ps_imm_rate = ps_imm_rate,
    n_island_species = n_island_species,
    n_mainland_species = n_mainland_species,
    carr_cap = carr_cap)

  expect_equal(expected, created)
})

test_that("use area constant diversity-dependent", {
  carr_cap <- 10
  ps_imm_rate <- 0.1
  n_island_species <- 5
  n_mainland_species <- 1
  hyper_pars <- DAISIE:::create_hyper_pars(0, 0)
  area <- 1
  created <- DAISIE:::get_immig_rate(
    gam = ps_imm_rate,
    A = area,
    num_spec = n_island_species,
    K = carr_cap,
    mainland_n = n_mainland_species)
  expected <- DAISIE_calc_clade_imm_rate(
    ps_imm_rate = ps_imm_rate,
    n_island_species = n_island_species,
    n_mainland_species = n_mainland_species,
    carr_cap = carr_cap)

  expect_equal(expected, created)
})

test_that("use area variable (ontogeny) diversity-dependent", {
  carr_cap <- 10
  ps_imm_rate <- 0.1
  n_island_species <- 5
  n_mainland_species <- 1
  hyper_pars <- DAISIE:::create_hyper_pars(0, 0)
  area <- 10
  created <- DAISIE:::get_immig_rate(
    gam = ps_imm_rate,
    A = area,
    num_spec = n_island_species,
    K = carr_cap,
    mainland_n = n_mainland_species)
  expected <- 0.095
  expect_equal(expected, created)
})
