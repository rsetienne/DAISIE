context("DAISIE_ML1")

test_that("DAISIE_ML1 works and simplex and subplex give the same answer", {
  skip_if(Sys.getenv("CI") == "", message = "Run only on CI")
  skip_on_cran()
  set.seed(42)
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
    methode = 'odeint::runge_kutta_cash_karp54',
    tol = c(1e-04, 1e-05, 1e-07),
    res = 15,
    tolint = c(1E-16, 1E-10),
    optimmethod = 'subplex',
    num_cycles = 4)
  expected_MLE <- data.frame(
    lambda_c = 1.8719078164871,
    mu = 1.93904173355875,
    K = Inf,
    gamma = 0.00742126301054098,
    lambda_a = 1.16344953861656,
    loglik = -76.7924216827578,
    df = 5L,
    conv = 0L
  )
  testthat::expect_equal(tested_MLE1, expected_MLE, tolerance = 1E-6)
  invisible(capture.output(tested_MLE2 <- DAISIE_ML1(
    datalist = datalist,
    initparsopt = as.numeric(tested_MLE1[1:5]),
    idparsopt = idparsopt,
    parsfix = parsfix,
    ddmodel = ddmodel,
    idparsfix = idparsfix,
    verbose = 0,
    methode = 'odeint::runge_kutta_cash_karp54',
    tol = c(1e-04, 1e-05, 1e-07),
    res = 15,
    tolint = c(1E-16, 1E-10),
    optimmethod = 'simplex',
    num_cycles = 1)))
  tested_MLE3 <- DAISIE_ML1(
    datalist = datalist,
    initparsopt = as.numeric(tested_MLE2[1:5]),
    idparsopt = idparsopt,
    parsfix = parsfix,
    ddmodel = ddmodel,
    idparsfix = idparsfix,
    verbose = 0,
    methode = 'odeint::runge_kutta_cash_karp54',
    tol = c(1e-04, 1e-05, 1e-07),
    res = 15,
    tolint = c(1E-16, 1E-10),
    optimmethod = 'subplex',
    num_cycles = 1)
  testthat::expect_equal(tested_MLE2, tested_MLE3, tolerance = 1E-6)
})

test_that("abuse", {
  #skip_if(Sys.getenv("CI") == "", message = "Run only on CI")
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
