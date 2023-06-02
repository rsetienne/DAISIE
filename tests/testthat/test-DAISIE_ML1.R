context("DAISIE_ML1")

test_that("use", {
  skip_if(Sys.getenv("CI") == "", message = "Run only on CI")
  skip_on_cran()
  data(Galapagos_datalist)
  datalist <- Galapagos_datalist
  initparsopt <- c(2.5, 2.7, 20, 0.009, 1.01)
  ddmodel <- 11
  idparsopt <- c(1,2,3,4,5)
  parsfix <- c()
  idparsfix <- c()
    tested_MLE <- DAISIE_ML1(
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
    lambda_c = 3.689104200780465,
    mu = 4.31030299415995,
    K = 906.6501180193454,
    gamma = 0.0173458887696076,
    lambda_a = 3.677789527566334,
    loglik = -64.2199684450019,
    df = 5L,
    conv = 0L
  )
  expect_equal(tested_MLE, expected_MLE)
})

test_that("abuse", {
  skip_if(Sys.getenv("CI") == "", message = "Run only on CI")
  skip_on_cran()
  utils::data(Galapagos_datalist)
  datalist <- Galapagos_datalist
  initparsopt <- c(2.5, 2.7, 20, 0.009, 1.01)
  ddmodel <- 11
  idparsopt <- 1:5
  parsfix <- NULL
  idparsfix <- NULL
  expect_error(
    DAISIE_ML1(
      datalist = "nonsense",
      initparsopt = initparsopt,
      idparsopt = idparsopt,
      parsfix = parsfix,
      ddmodel = ddmodel,
      idparsfix = idparsfix,
      verbose = 0
    )
  )
})
