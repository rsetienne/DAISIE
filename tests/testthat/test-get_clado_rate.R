context("get_clado_rate")

test_that("use area constant diversity-independent", {
  ps_clado_rate <- 0.2
  carr_cap <- Inf
  n_species <- 4
  hyper_pars <- DAISIE:::create_hyper_pars(d = 0, x = 0)
  area <- 1
  created <- DAISIE:::get_clado_rate(
    lac = ps_clado_rate,
    hyper_pars = hyper_pars,
    num_spec = n_species,
    K = carr_cap,
    A = area
  )
  expected <- DAISIE_calc_clade_clado_rate(
    ps_clado_rate = ps_clado_rate,
    n_species = n_species,
    carr_cap = carr_cap
  )
  expect_equal(created, expected)
})

test_that("use area constant diversity-dependent", {
  ps_clado_rate <- 0.2
  carr_cap <- 9
  n_species <- 4
  hyper_pars <- DAISIE:::create_hyper_pars(d = 0, x = 0)
  area <- 1
  created <- DAISIE:::get_clado_rate(
    lac = ps_clado_rate,
    hyper_pars = hyper_pars,
    num_spec = n_species,
    K = carr_cap,
    A = area
  )
  expected <- DAISIE_calc_clade_clado_rate(
    ps_clado_rate = ps_clado_rate,
    n_species = n_species,
    carr_cap = carr_cap
  )
  expect_equal(created, expected)
})

test_that("use area variable diversity-independent", {
  ps_clado_rate <- 0.2
  carr_cap <- Inf
  n_species <- 4
  hyper_pars <- DAISIE:::create_hyper_pars(d = 0.2, x = 0.1)
  area <- 10
  created <- DAISIE:::get_clado_rate(
    lac = ps_clado_rate,
    hyper_pars = hyper_pars,
    num_spec = n_species,
    K = carr_cap,
    A = area
  )
  expected <- 1.267914553968891
  expect_equal(created, expected)
})

test_that("use area variable diversity-dependent", {
  ps_clado_rate <- 0.2
  carr_cap <- 9
  n_species <- 4
  hyper_pars <- DAISIE:::create_hyper_pars(d = 0.2, x = 0.1)
  area <- 10
  created <- DAISIE:::get_clado_rate(
    lac = ps_clado_rate,
    hyper_pars = hyper_pars,
    num_spec = n_species,
    K = carr_cap,
    A = area
  )
  expected <- 1.211562796014718
  expect_equal(created, expected)
})
