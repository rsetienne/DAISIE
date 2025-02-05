test_that("DAISIE_logp0 is correct", {
  pars1 <- c(2.000, 2.700, 20.000, 0.009, 1.010)
  pars2 <- c(1.0e+02, 1.1e+01, 0.0e+00, 0.0e+00, NA, 0.0e+00, 1.0e-04,
             1.0e-05, 1.0e-07, 3.0e+03, 9.5e-01, 9.8e-01)
  island_age <- 4
  logp0 <- DAISIE_logp0(pars1 = pars1,
                        pars2 = pars2,
                        island_age = island_age)
  testthat::expect_true(is.numeric(logp0))
  testthat::expect_equal(logp0, -0.005413629)
})
