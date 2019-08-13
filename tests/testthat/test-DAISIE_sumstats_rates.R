context("DAISIE_sumstats_rates")

test_that("use", {
  out <- DAISIE_calc_sumstats_pcrates()
  expect_true(is.list(out))
  
})

test_that("abuse", {
  expect_error()
})
