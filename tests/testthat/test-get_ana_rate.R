test_that("use", {
  ps_ana_rate <- 1
  n_immigrants <- 5
  ana_rate <- get_ana_rate(
    laa = ps_ana_rate,
    num_immigrants = n_immigrants)
  created <- get_ana_rate(
    laa = 1,
    num_immigrants = 5)
  expected <- ps_ana_rate * n_immigrants
  testthat::expect_equal(expected, created)
})
