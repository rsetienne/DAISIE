test_that("use", {
  daisie_params <- create_test_daisie_params()
  expect_equal(daisie_params$time, 3)
  expect_equal(daisie_params$M, 1)
  expect_equal(daisie_params$pars, c(2.5, 2.6, Inf, 0.01, 1.0))
  expect_equal(daisie_params$replicates, 1)
})
