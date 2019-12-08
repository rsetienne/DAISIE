context("get_ana_rate")

test_that("ana_rate is a number", {
  expect_silent(
    is.numeric(
      get_ana_rate(
        laa = 1,
        hyper_pars = c(1, 1, 1, 1),
        dist_pars = 1,
        num_immigrants = 10
      )
    )
  )
})

test_that("classic behaviour", {
  ps_ana_rate <- 1
  n_immigrants <- 5
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


  test_that("use hyper_pars", {
    expect_silent(
      is.numeric(
        get_ana_rate(
          laa = 1,
          hyper_pars = c(1, 1, 1, 1),
          dist_pars = 1,
          num_immigrants = 5
        )
      )
    )
  })
