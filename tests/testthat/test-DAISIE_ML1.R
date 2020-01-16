context("test-DAISIE_ML1")

test_that("use", {
  if (Sys.getenv("TRAVIS") != "") {

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
      verbose = 0
    )
    expected_MLE <- data.frame(
      lambda_c = 2.55847849219339,
      mu = 2.68768191590176,
      K = 6765.0637400135,
      gamma = 0.00932987953669849,
      lambda_a = 1.00838182578826,
      loglik = -76.0001379108545,
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
