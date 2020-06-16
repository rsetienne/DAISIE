context("DAISIE_ML1")

test_that("use", {
  if (Sys.getenv("TRAVIS") != "") {

    # This is a rough MLE test, built for fast execution. A more thorough test
    # can be found in the GitHub repository Neves-P/DAISIEtesting

    utils::data(Galapagos_datalist)
    datalist <- Galapagos_datalist
    initparsopt <- c(2.5, 2.7, 20, 0.009, 1.01)
    ddmodel <- 11
    idparsopt <- 1:5
    parsfix <- NULL
    idparsfix <- NULL
    tested_MLE <- DAISIE:::DAISIE_ML1(
      datalist = datalist,
      initparsopt = initparsopt,
      idparsopt = idparsopt,
      parsfix = parsfix,
      ddmodel = ddmodel,
      idparsfix = idparsfix,
      verbose = 0,
      tol = c(0.01, 0.1, 0.001),
      res = 15,
      tolint = c(0.1, 0.01)
    )
    expected_MLE <- data.frame(
      lambda_c = 4.1558928147225656,
      mu = 4.9989880317542603,
      K = 514.14066875326353,
      gamma = 0.020398370167867434,
      lambda_a = 3.7589983231693607,
      loglik = -63.993218343764838,
      df = 5L,
      conv = 0L
    )
    expect_equal(tested_MLE, expected_MLE)
  } else {
    skip("Run only on Travis")
  }

})

test_that("abuse", {
  if (Sys.getenv("TRAVIS") != "") {
    utils::data(Galapagos_datalist)
    datalist <- Galapagos_datalist
    initparsopt <- c(2.5, 2.7, 20, 0.009, 1.01)
    ddmodel <- 11
    idparsopt <- 1:5
    parsfix <- NULL
    idparsfix <- NULL
    expect_error(
      DAISIE:::DAISIE_ML1(
        datalist = "nonsense",
        initparsopt = initparsopt,
        idparsopt = idparsopt,
        parsfix = parsfix,
        ddmodel = ddmodel,
        idparsfix = idparsfix,
        verbose = 0
      )
    )
  } else {
    skip("Run only on Travis")
  }
})
