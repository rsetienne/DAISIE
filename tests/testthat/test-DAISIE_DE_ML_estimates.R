

test_that("DAISIE_ML_CS: DAISIE_DE with equal_extinction = TRUE matches DAISIE", {
  #skip("WIP")
  utils::data(Galapagos_datalist)

  ML_estimates_DAISIE <- DAISIE_ML_CS(
    datalist = Galapagos_datalist,
    initparsopt = c(2.280147, 2.421845, 0.00391762, 1.313127),
    idparsopt = c(1, 2, 4, 5),
    parsfix = Inf,
    idparsfix = 3,
    ddmodel = 0,
    function_to_optimize = "DAISIE"
  )

  ML_estimates_DAISIE_DE <- DAISIE_ML_CS(
    datalist = Galapagos_datalist,
    initparsopt = c(2.280147, 2.421845, 0.00391762, 1.313127),
    idparsopt = c(1, 2, 4, 5),
    parsfix = Inf,
    idparsfix = 3,
    ddmodel = 0,
    CS_version = list(model = 1, function_to_optimize = 'DAISIE_DE'),
    equal_extinction = TRUE
  )

  testthat::expect_equal(ML_estimates_DAISIE_DE$loglik, ML_estimates_DAISIE$loglik)
})
