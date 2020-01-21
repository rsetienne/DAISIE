context("get_ana_rate")

test_that("use area constant diversity-independent without hyper_pars", {
  ps_ana_rate <- 1
  n_immigrants <- 5
  default_pars <- create_default_pars(
      island_ontogeny = 0,
      sea_level = 0,
      totaltime = 1,
      hyper_pars = NULL,
      dist_pars = NULL,
      ext_pars = NULL,
      pars = c(0, 0, 0, 0, 0),
      area_pars = NULL
    )
  expect_silent(ana_rate <- DAISIE:::get_ana_rate(
    laa = 1,
    hyper_pars = default_pars$hyper_pars,
    dist_pars = default_pars$dist_pars,
    num_immigrants = 5)
  )
  expect_true(is.numeric(ana_rate))
  expected <- DAISIE_calc_clade_ana_rate(
    ps_ana_rate = ps_ana_rate,
    n_immigrants = n_immigrants
  )
  created <- get_ana_rate(
    laa = 1,
    hyper_pars = default_pars$hyper_pars,
    dist_pars = default_pars$dist_pars,
    num_immigrants = 5
  )
  expect_equal(expected, created)
})

test_that("use area constant diversity-independent with hyper_pars", {
  ps_ana_rate <- 1
  n_immigrants <- 5
  hyper_pars <- create_hyper_pars(d_0 = 1, x =  1, alpha =  1, beta =  1)
  dist_pars <- create_dist_pars(D = 1)
  expect_silent(ana_rate <- DAISIE:::get_ana_rate(
    laa = ps_ana_rate,
    hyper_pars = hyper_pars,
    dist_pars = dist_pars,
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
    hyper_pars = hyper_pars,
    dist_pars = dist_pars,
    num_immigrants = 5
  )
  expect_equal(expected, created)
})
