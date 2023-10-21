test_that("use area constant", {
  ps_ext_rate <- 2
  n_species <- 4
  hyper_pars = create_hyper_pars(d = 0, x = 0)
  area <- 1
  created <- get_ext_rate(
    mu = ps_ext_rate,
    hyper_pars = hyper_pars,
    extcutoff = 1000,
    num_spec = n_species,
    A = area)
  expected <- ps_ext_rate * n_species

  testthat::expect_equal(created, expected)
})

test_that("use area variable", {
  ps_ext_rate <- 2
  n_species <- 4
  hyper_pars = create_hyper_pars(d = 0.2, x = 0.1)
  area <- 10
  created <- get_ext_rate(
    mu = ps_ext_rate,
    hyper_pars = hyper_pars,
    extcutoff = 1000,
    num_spec = n_species,
    A = area)
  expected <- 6.354625877794252
  testthat::expect_equal(created, expected)
})
