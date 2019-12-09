context("get_ana_rate")

test_that("use area constant diversity-independent without hyper_pars", {
  ps_ana_rate <- 1
  n_immigrants <- 5
  expect_silent(ana_rate <- DAISIE:::get_ana_rate(
    laa = 1,
    hyper_pars = NULL,
    dist_pars = NULL,
    num_immigrants = 5)
  )
  expect_true(is.numeric(ana_rate))
  expected <- DAISIE_calc_clade_ana_rate(
    ps_ana_rate = ps_ana_rate,
    n_immigrants = n_immigrants
  )
  created <- get_ana_rate(
    laa = 1,
    hyper_pars = NULL,
    dist_pars = NULL,
    num_immigrants = 5
  )
  expect_equal(expected, created)
})

test_that("use area constant diversity-independent with hyper_pars", {
  ps_ana_rate <- 1
  n_immigrants <- 5
  expect_silent(ana_rate <- DAISIE:::get_ana_rate(
    laa = 1,
    hyper_pars = c(1, 1, 1, 1),
    dist_pars = 1,
    num_immigrants = 5
  )
  )
  expect_true(is.numeric(ana_rate))
  expected <- DAISIE_calc_clade_ana_rate(
    ps_ana_rate = ps_ana_rate,
    n_immigrants = n_immigrants
  )
  created <- get_ana_rate(
    laa = 1,
    hyper_pars = c(1,1,1,1),
    dist_pars = 1,
    num_immigrants = 5
  )
  expect_equal(expected, created)
})
