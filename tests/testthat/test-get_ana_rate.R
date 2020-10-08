context("get_ana_rate")

test_that("use", {
  ps_ana_rate <- 1
  n_immigrants <- 5
  ana_rate <- DAISIE:::get_ana_rate(
    laa = ps_ana_rate,
    num_immigrants = n_immigrants)
  created <- DAISIE:::get_ana_rate(
    laa = 1,
    num_immigrants = 5)
  expected <- DAISIE_calc_clade_ana_rate(
    ps_ana_rate = ps_ana_rate,
    n_immigrants = n_immigrants)
  expect_equal(expected, created)
})
