test_that("DAISIE_ML_CS: DAISIE_DE with equal_extinction = TRUE matches DAISIE", {
  #skip("WIP")
  utils::data(Galapagos_datalist)

  ML_estimates_DAISIE <- DAISIE_ML_CS(
    datalist = Galapagos_datalist,
    initparsopt = c(2.550682, 2.683817, 0.009344, 1.00728),
    idparsopt = c(1, 2, 4, 5),
    parsfix = Inf,
    idparsfix = 3,
    ddmodel = 0,
    verbose = 0,
    CS_version = list(model = 1,function_to_optimize = "DAISIE")
  )

  ML_estimates_DAISIE_DE <- DAISIE_ML_CS(
    datalist = Galapagos_datalist,
    initparsopt = c(2.550682, 2.683817, 0.009344, 1.00728),
    idparsopt = c(1, 2, 4, 5),
    parsfix = Inf,
    idparsfix = 3,
    ddmodel = 0,
    verbose = 0,
    CS_version = list(model = 1, function_to_optimize = 'DAISIE_DE'),
    equal_extinction = TRUE
  )

  testthat::expect_equal(ML_estimates_DAISIE_DE$loglik, ML_estimates_DAISIE$loglik, tol = 1E-6)
  testthat::expect_equal(ML_estimates_DAISIE_DE, ML_estimates_DAISIE, tol = 1E-3)
})
