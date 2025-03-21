test_that("use", {
  skip_if(Sys.getenv("CI") == "", message = "Run only on CI")
  skip_on_cran()
  utils::data(Macaronesia_datalist, package = "DAISIE")
    tested_MLE <- DAISIE_ML2(
      datalist = Macaronesia_datalist,
      initparsopt = c(
        1.053151832,
        0.052148979,
        0.512939011,
        0.133766934,
        0.152763179
      ),
      idparsmat = rbind(
        1:5,
        c(6, 2, 3, 7, 5),
        1:5,1:5
      ),
      idparsopt = c(2, 4, 5, 6, 7),
      parsfix = c(0, Inf),
      idparsfix = c(1, 3),
      methode = 'odeint::runge_kutta_cash_karp54',
      tol = c(0.01, 0.1, 0.001),
      res = 15,
      tolint = c(1E-16, 1E-10),
      verbose = 0
    )
  #expected_MLE <- data.frame(
  #  lambda_c = c(0.0,
  #               0.1198600880941345,
  #               0.0,
  #               0.0),
  #  mu = rep(1.0604828698048443,4),
  #  K = rep(Inf,4),
  #  gamma = c(0.0544106971300572,
  #            0.1706160417441286,
  #            0.0544106971300572,
  #            0.0544106971300572),
  #  lambda_a = rep(0.4749396123469573,4),
  #  loglik = rep(-449.6525311168169310,4),
  #  df = rep(5L,4),
  #  conv = rep(0L,4)
  #)
  expected_MLE <- data.frame(
    lambda_c = c(0,0.1168013973697136,0,0),
    mu = rep(1.063126398848375,4),
    K = rep(Inf,4),
    gamma = c(0.05379104379285973,0.17160055800328039,0.05379104379285973,0.05379104379285973),
    lambda_a = 0.4675764819674205,
    loglik = rep(-449.9125224959051,4),
    df = rep(5L,4),
    conv = rep(0L,4)
  )
  testthat::expect_equal(tested_MLE, expected_MLE)
})


test_that("abuse", {
  skip_if(Sys.getenv("CI") == "", message = "Run only on CI")
  skip_on_cran()
  testthat::expect_error(tested_MLE <- DAISIE_ML2(
    datalist = "nonsense",
    initparsopt = c(
      1.053151832,
      0.052148979,
      0.512939011,
      0.133766934,
      0.152763179
    ),
    idparsmat = rbind(
      1:5,
      c(6, 2, 3, 7, 5),
      1:5,1:5
    ),
    idparsopt = c(2, 4, 5, 6, 7),
    parsfix = c(0, Inf),
    idparsfix = c(1, 3),
    tol = c(0.01, 0.1, 0.001),
    res = 15,
    tolint = c(0.1, 0.01)
  ))
})
