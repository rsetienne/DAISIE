context("DAISIE_loglik_CS_choice")

test_that("DAISIE_loglik_CS_choice produces correct output for CS_version 1", {
  pars1 <- c(2.000, 2.700, 20.000, 0.009, 1.010)
  pars2 <- c(1.0e+02, 1.1e+01, 0.0e+00, 0.0e+00, NA, 0.0e+00, 1.0e-04,
             1.0e-05, 1.0e-07, 3.0e+03, 9.5e-01, 9.8e-01)
  brts <- c(4.0000, 3.0282, 1.3227, 0.8223, 0.4286, 0.3462, 0.2450, 0.0808,
            0.0527, 0.0327, 0.0221, 0.1180, 0.0756, 0.0525, 0.0322, 0.0118)
  stac <- 2
  missnumspec <- 0
  loglik <- DAISIE_loglik_CS_choice(pars1 = pars1,
                                    pars2 = pars2,
                                    brts = brts,
                                    stac = stac,
                                    missnumspec = missnumspec)
  expect_true(is.numeric(loglik))
  expect_equal(loglik, -17.6550433826)
})

test_that("DAISIE_loglik_CS_choice produces correct output for CS_version 2", {
  pars1 <- c(2.000, 2.700, 20.000, 0.009, 1.010)
  pars2 <- c(1.0e+02, 1.1e+01, 0.0e+00, 0.0e+00, NA, 0.0e+00, 1.0e-04,
             1.0e-05, 1.0e-07, 3.0e+03, 9.5e-01, 9.8e-01)
  brts <- c(4.0000, 3.0282, 1.3227, 0.8223, 0.4286, 0.3462, 0.2450, 0.0808,
            0.0527, 0.0327, 0.0221, 0.1180, 0.0756, 0.0525, 0.0322, 0.0118)
  stac <- 2
  missnumspec <- 0
  CS_version <- create_CS_version(model = 2,
                                  pick_parameter = "cladogenesis",
                                  distribution = "gamma",
                                  sd = 1)
  loglik <- DAISIE_loglik_CS_choice(pars1 = pars1,
                                    pars2 = pars2,
                                    brts = brts,
                                    stac = stac,
                                    missnumspec = missnumspec,
                                    CS_version = CS_version)
  expect_true(is.numeric(loglik))
  expect_equal(loglik, -9.55117524011)
})

test_that("DAISIE_loglik_CS_choice produces correct output for CS_version 0", {
  pars1 <- c(2.000, 2.700, 20.000, 0.009, 1.010)
  pars2 <- c(1.0e+02, 1.1e+01, 0.0e+00, 0.0e+00, NA, 0.0e+00, 1.0e-04,
             1.0e-05, 1.0e-07, 3.0e+03, 9.5e-01, 9.8e-01)
  brts <- c(4.0000, 3.0282, 1.3227, 0.8223, 0.4286, 0.3462, 0.2450, 0.0808,
            0.0527, 0.0327, 0.0221, 0.1180, 0.0756, 0.0525, 0.0322, 0.0118)
  stac <- 2
  missnumspec <- 0
  CS_version <- 0
  loglik <- DAISIE_loglik_CS_choice(pars1 = pars1,
                                    pars2 = pars2,
                                    brts = brts,
                                    stac = stac,
                                    missnumspec = missnumspec,
                                    CS_version = CS_version)
  expect_true(is.numeric(loglik))
  expect_equal(loglik, -17.5608831694)
})
