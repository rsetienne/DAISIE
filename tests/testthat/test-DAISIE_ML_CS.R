test_that("relaxed-rate DAISIE_ML_CS produces correct output", {
  skip("Too slow to run")
  utils::data(Galapagos_datalist)
  CS_version <- create_CS_version(model = 2,
                                  relaxed_par = "cladogenesis")
  RR_clado <- DAISIE_ML_CS(
    datalist = Galapagos_datalist,
    initparsopt = c(2.5, 2.7, 20, 0.009, 1.01, 2),
    idparsopt = 1:6,
    parsfix = NULL,
    idparsfix = NULL,
    ddmodel = 11,
    CS_version = CS_version,
    methode = 'odeint::runge_kutta_cash_karp54',
    verbose = 1)
  testthat::expect_true(is.numeric(RR_clado))
  testthat::expect_true(is.numeric(result$loglik))
})

test_that("relaxed-rate DAISIE_ML_CS produces correct output using simplex", {
  skip("Too slow to run")
  utils::data(Galapagos_datalist)
  CS_version <- create_CS_version(model = 2,
                                  relaxed_par = "cladogenesis")
  RR_clado <- DAISIE_ML_CS(
    datalist = Galapagos_datalist,
    initparsopt = c(2.5, 2.7, 20, 0.009, 1.01, 2),
    idparsopt = 1:6,
    parsfix = NULL,
    idparsfix = NULL,
    ddmodel = 11,
    CS_version = CS_version,
    optimmethod = "simplex")
  testthat::expect_true(is.numeric(RR_clado))
  testthat::expect_true(is.numeric(result$loglik))
})

test_that("multi-rate DAISIE_ML_CS converges to constant rate", {
  skip("WIP")
  utils::data(Galapagos_datalist)
  CR <- DAISIE_ML_CS(datalist = Galapagos_datalist,
                     initparsopt = c(2.5, 2.7, 20, 0.009, 1.01),
                     idparsopt = 1:5,
                     parsfix = NULL,
                     idparsfix = NULL,
                     ddmodel = 11,
                     CS_version = list(model = 1, function_to_optimize = 'DAISIE'))

  utils::data(Galapagos_datalist)
  CS_version <- create_CS_version(model = 2,
                                  relaxed_par = "cladogenesis")
  RR_clado <- DAISIE_ML_CS(datalist = Galapagos_datalist,
                                initparsopt = c(2.5, 2.7, 20, 0.009, 1.01),
                                idparsopt = 1:5,
                                parsfix = NULL,
                                idparsfix = NULL,
                                ddmodel = 11,
                                CS_version = CS_version)
  testthat::expected_equal(CR, RR_clado)
})

