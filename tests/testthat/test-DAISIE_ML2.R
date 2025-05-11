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
  expected_MLE <- data.frame(
    lambda_c = c(0,0.116801397369714,0,0),
    mu = rep(1.06312639884837,4),
    K = rep(Inf,4),
    gamma = c(0.0537910437928597,0.1716005580032804,0.0537910437928597,0.0537910437928597),
    lambda_a = rep(0.46757648196742,4),
    loglik = rep(-449.913214736693,4),
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
