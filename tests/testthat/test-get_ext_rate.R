context("get_ext_rate")

test_that("use area constant", {
  ps_ext_rate <- 2
  n_species <- 4
  hyper_pars = DAISIE:::create_hyper_pars(d = 0, x = 0)
  area <- 1
  created <- DAISIE:::get_ext_rate(
    mu = ps_ext_rate,
    hyper_pars = hyper_pars,
    extcutoff = 1000,
    num_spec = n_species,
    A = area)
  expected <- DAISIE_calc_clade_ext_rate(
    ps_ext_rate = ps_ext_rate,
    n_species = n_species)
  expect_equal(created, expected)
})

test_that("use area variable", {
  ps_ext_rate <- 2
  n_species <- 4
  hyper_pars = DAISIE:::create_hyper_pars(d = 0.2, x = 0.1)
  area <- 10
  created <- DAISIE:::get_ext_rate(
    mu = ps_ext_rate,
    hyper_pars = hyper_pars,
    extcutoff = 1000,
    num_spec = n_species,
    A = area)
  expected <- 6.354625877794252
  expect_equal(created, expected)
})
