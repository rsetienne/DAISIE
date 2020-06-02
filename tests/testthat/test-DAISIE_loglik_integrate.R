context("DAISIE_loglik_integrate")

test_that("DAISIE_loglik_integrate produces correct ouput on single lineage", {
  pars1 <- c(2.000, 2.700, 20.000, 0.009, 1.010)
  pars2 <- c(1.0e+02, 1.1e+01, 0.0e+00, 0.0e+00, NA, 0.0e+00, 1.0e-04,
             1.0e-05, 1.0e-07, 3.0e+03, 9.5e-01, 9.8e-01)
  brts <- 4
  stac <- 0
  missnumspec <- 0
  methode <- "lsodes"
  CS_version <- create_CS_version(model = 2,
                                  pick_parameter = 'carrying_capacity',
                                  distribution = 'gamma',
                                  sd = 2,
                                  multi_rate_optim_method = 'optimize' )
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
  testthat::expect_equal(loglik, -0.005410626, tolerance = 1E-5)

  pars1 <- c(0.25, 0.1, 300, 0.05, 1)
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
  testthat::expect_equal(loglik, -38.5417852385957644, tolerance = 1E-5)

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
  CS_version <- create_CS_version(model = 2,
                                  pick_parameter = 'carrying_capacity',
                                  distribution = 'gamma',
                                  sd = 2,
                                  multi_rate_optim_method = 'optimize')
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
  testthat::expect_equal(loglik, -17.2793971737510610, tolerance = 1E-5)

  pars1 <- c(0.25, 0.1, 300, 0.05, 1)
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
  testthat::expect_equal(loglik, -46.1179137501893877, tolerance = 1E-5)

})
