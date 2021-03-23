context("DAISIE_loglik_integrate")

test_that("DAISIE_loglik_integrate produces correct ouput on single lineage", {
  pars1 <- c(2.000, 2.700, 20.000, 0.009, 1.010)
  pars2 <- c(1.0e+02, 1.1e+01, 0.0e+00, 0.0e+00, NA, 0.0e+00, 1.0e-04,
             1.0e-05, 1.0e-07, 3.0e+03, 9.5e-01, 9.8e-01)
  brts <- 4
  stac <- 0
  missnumspec <- 0
  methode <- "lsodes"
  CS_version <- list(model = 2,
                     relaxed_par = 'carrying_capacity',
                     sd = 2)
  abstolint <- 1e-16
  reltolint <- 1e-10
  verbose <- FALSE
  loglik <- DAISIE_loglik_integrate(
    pars1 = pars1,
    pars2 = pars2,
    brts = brts,
    stac = stac,
    missnumspec = missnumspec,
    methode = methode,
    CS_version = CS_version,
    abstolint = abstolint,
    reltolint = reltolint,
    verbose = verbose
  )
  testthat::expect_equal(loglik, -0.00541062585063765)
})

test_that("DAISIE_loglik_integrate produces correct ouput on radiation", {
  pars1 <- c(2.000, 2.700, 20.000, 0.009, 1.010)
  pars2 <- c(1.0e+02, 1.1e+01, 0.0e+00, 0.0e+00, NA, 0.0e+00, 1.0e-04,
             1.0e-05, 1.0e-07, 3.0e+03, 9.5e-01, 9.8e-01)
  brts <- c(4.0000, 3.0282, 1.3227, 0.8223, 0.4286, 0.3462, 0.2450, 0.0808,
            0.0527, 0.0327, 0.0221, 0.1180, 0.0756, 0.0525, 0.0322, 0.0118)
  stac <- 2
  missnumspec <- 0
  methode <- "lsodes"
  CS_version <- list(model = 2,
                     relaxed_par = 'carrying_capacity',
                     sd = 10)
  abstolint <- 1e-16
  reltolint <- 1e-10
  verbose <- FALSE
  loglik <- DAISIE_loglik_integrate(
    pars1 = pars1,
    pars2 = pars2,
    brts = brts,
    stac = stac,
    missnumspec = missnumspec,
    methode = methode,
    CS_version = CS_version,
    abstolint = abstolint,
    reltolint = reltolint,
    verbose = verbose
  )
  testthat::expect_equal(loglik, -15.1289048939324)
})

test_that("DAISIE_loglik_integrand produces correct output", {
  output <- DAISIE_loglik_integrand(
    DAISIE_par = 1,
    pars1 = c(2.550687345, 2.683454548, 10.000000000,
              0.009332070, 1.010073119),
    pars2 = c(100, 0, 0, 0, NA),
    brts = 4,
    stac = 0,
    missnumspec = 0,
    methode = "lsodes",
    abstolint = 1e-16,
    reltolint = 1e-10,
    verbose = FALSE,
    pick = 1,
    par_mean = 2.550687345,
    par_sd = 1)
  expect_equal(output, -2.13638048160996)
})

test_that("rho produces correct output", {
  output <- rho(DAISIE_par = 0.5,
                DAISIE_dist_pars = list(
                  par_mean = 1,
                  par_sd = 1))
  expect_equal(output, -0.5)
})

test_that("transform_gamma_pars produces correct output", {
  output <- transform_gamma_pars(par_mean = 1, par_sd = 1)
  expect_equal(output, list(shape = 1, scale = 1))
})
