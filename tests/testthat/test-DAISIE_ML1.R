context("DAISIE_ML1")

test_that("DAISIE_ML1 works and simplex and subplex give the same answer", {
  skip_if(Sys.getenv("CI") == "", message = "Run only on CI")
  skip_on_cran()
  data(Galapagos_datalist)
  datalist <- Galapagos_datalist
  initparsopt <- c(2.5, 2.7, 20, 0.009, 1.01)
  ddmodel <- 11
  idparsopt <- c(1,2,3,4,5)
  parsfix <- c()
  idparsfix <- c()
    tested_MLE1 <- DAISIE_ML1(
      datalist = datalist,
      initparsopt = initparsopt,
      idparsopt = idparsopt,
      parsfix = parsfix,
      ddmodel = ddmodel,
      idparsfix = idparsfix,
      verbose = 0,
      tol = c(0.0001, 0.001, 0.0001),
      res = 15,
      tolint = c(0.1, 0.01),
      optimmethod = 'subplex',
      num_cycles = 3)
    tested_MLE2 <- DAISIE_ML1(
      datalist = datalist,
      initparsopt = as.numeric(tested_MLE1[1:5]),
      idparsopt = idparsopt,
      parsfix = parsfix,
      ddmodel = ddmodel,
      idparsfix = idparsfix,
      verbose = 0,
      tol = c(0.0001, 0.001, 0.0001),
      res = 15,
      tolint = c(0.1, 0.01),
      optimmethod = 'simplex')
  testthat::expect_equal(tested_MLE1, tested_MLE2)
  expected_MLE <- data.frame(
    lambda_c = 4.0275356252375420,
    mu = 4.8740259531255852,
    K = 2001.9278290878507960,
    gamma = 0.0202054262850127,
    lambda_a = 4.0124628182465916,
    loglik = -63.9341963365704231,
    df = 5L,
    conv = 0L
  )
  testthat::expect_equal(tested_MLE1, expected_MLE)
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
  testthat::expect_error(
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
