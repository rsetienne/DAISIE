test_that("DAISIE_ML_CS: DAISIE_DE with equal_extinction = TRUE matches DAISIE", {
  #skip("WIP")
  utils::data(Galapagos_datalist)

  invisible(capture.output(ML_estimates_DAISIE <- DAISIE_ML_CS(
    datalist = Galapagos_datalist,
    initparsopt = c(2.550682, 2.683817, 0.009344, 1.00728),
    idparsopt = c(1, 2, 4, 5),
    parsfix = Inf,
    idparsfix = 3,
    ddmodel = 0,
    verbose = 0,
    CS_version = list(
      model = 1,
      function_to_optimize = "DAISIE")
  )))

  invisible(capture.output(ML_estimates_DAISIE_DE <- DAISIE_ML_CS(
    datalist = Galapagos_datalist,
    initparsopt = c(2.550682, 2.683817, 0.009344, 1.00728),
    idparsopt = c(1, 2, 4, 5),
    parsfix = Inf,
    idparsfix = 3,
    ddmodel = 0,
    verbose = 0,
    methode = 'lsodes',
    CS_version = list(model = 1, function_to_optimize = 'DAISIE_DE'),
    equal_extinction = TRUE
  )))

  testthat::expect_equal(ML_estimates_DAISIE_DE$loglik, ML_estimates_DAISIE$loglik, tol = 1E-6)
  testthat::expect_equal(ML_estimates_DAISIE_DE, ML_estimates_DAISIE, tol = 1E-3)
})

test_that("DAISIE_DE and DAISIE give same results when there are missing species", {

  pars1 <- c(0.2, 0.1, 0.1, 0.02, 0.03)
  brts <- c(4.000, 0.855)
  missnumspec <- 5
  loglik_DE <- DAISIE_DE_logpES(brts = brts,
                                missnumspec = missnumspec,
                                pars1 = pars1,
                                methode = "lsodes",
                                reltolint = 1e-16,
                                abstolint = 1e-16)
  pars1[3] <- Inf
  lik <- 0
  S <- length(brts) - 1
  for(i in 0:300) {
    lik <- lik + exp(dbinom(i, S + i,prob = missnumspec/(S + missnumspec), log = TRUE) +
                       DAISIE:::DAISIE_loglik(brts = brts,
                                              pars2 = c(150, 0, 0, 0),
                                              pars1 = pars1,
                                              stac = 2,
                                              missnumspec = i,
                                              methode = "odeint::runge_kutta_fehlberg78",
                                              reltolint = 1e-16,
                                              abstolint = 1e-16))
  }
  loglik <- log(lik)
  testthat::expect_equal(loglik, loglik_DE)
})
