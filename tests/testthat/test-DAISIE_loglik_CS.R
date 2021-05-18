context("DAISIE_loglik_CS")

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
  expect_equal(loglik, -17.6535269346579)
})

test_that("DAISIE_loglik_CS_choice produces correct output for relaxed-rate
          model (CS_version = 2)", {
            skip("Produces DLSODES warnings")
  pars1 <- c(2.000, 2.700, 20.000, 0.009, 1.010)
  pars2 <- c(1.0e+02, 1.1e+01, 0.0e+00, 0.0e+00, NA, 0.0e+00, 1.0e-04,
             1.0e-05, 1.0e-07, 3.0e+03, 9.5e-01, 9.8e-01)
  brts <- c(4.0000, 3.0282, 1.3227, 0.8223, 0.4286, 0.3462, 0.2450, 0.0808,
            0.0527, 0.0327, 0.0221, 0.1180, 0.0756, 0.0525, 0.0322, 0.0118)
  stac <- 2
  missnumspec <- 0
  CS_version <- list(model = 2,
                     relaxed_par = "cladogenesis",
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

test_that("DAISIE_loglik_all produces correct output for relaxed-rate model", {
  skip("Produces DLSODES warnings")
  utils::data(Galapagos_datalist)
  loglik <- DAISIE::DAISIE_loglik_all(
    pars1 = c(2.55068735, 2.68345455, 10.00000000, 0.00933207, 1.01007312),
    pars2 = c(100, 0, 0, 0, NA),
    datalist = Galapagos_datalist,
    methode = "lsodes",
    CS_version = list(model = 2,
                      relaxed_par = "cladogenesis",
                      sd = 1),
    abstolint = 1e-16,
    reltolint = 1e-10
  )
  expect_true(is.numeric(loglik))
  expect_equal(loglik, --77.5108137039949)
})

test_that("DAISIE_loglik produces correct output", {
  output <- DAISIE_loglik(pars1 = c(2.061154e-09, 2.683455e+00, 1.000000e+01,
                                    9.332070e-03, 1.010073e+00),
                            pars2 = c(100, 0, 0, 0, NA),
                            brts = 4,
                            stac = 0,
                            missnumspec = 0,
                            methode = "lsodes",
                            abstolint = 1E-16,
                            reltolint = 1E-10,
                            verbose = FALSE)
  expect_equal(output, -0.00347317077256095)
})
