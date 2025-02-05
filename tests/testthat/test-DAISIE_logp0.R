test_that("DAISIE_logp0 is correct", {
  pars1 <- c(2.000, 2.700, 20.000, 0.009, 1.010)
  pars2 <- c(1.0e+02, 1.1e+01)
  island_age <- 4
  logp0 <- DAISIE_logp0(pars1 = pars1,
                        pars2 = pars2,
                        island_age = island_age)
  testthat::expect_true(is.numeric(logp0))
  testthat::expect_equal(logp0, -0.005413629)
})
