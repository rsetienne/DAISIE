context("DAISIE_calc_clade_ana_rate")

test_that("use", {
  ps_ana_rate <- 1.2
  n_immigrants <- 42
  created <- DAISIE_calc_clade_ana_rate(
    ps_ana_rate = ps_ana_rate,
    n_immigrants = n_immigrants
  )
  expected <- ps_ana_rate * n_immigrants
  expect_equal(created, expected)
})

test_that("abuse", {
  ps_ana_rate <- 1.2
  n_immigrants <- 42
  expect_error(
    DAISIE_calc_clade_ana_rate(
      ps_ana_rate = -123.456,
      n_immigrants = n_immigrants
    )
  )
  expect_error(
    DAISIE_calc_clade_ana_rate(
      ps_ana_rate = ps_ana_rate,
      n_immigrants = -123456
    )
  )
})
