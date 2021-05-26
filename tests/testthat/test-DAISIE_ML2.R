context("DAISIE_ML2")
test_that("use", {

  skip_if(Sys.getenv("CI") == "", message = "Run only on CI")
  utils::data(Macaronesia_datalist, package = "DAISIE")
  invisible(capture.output(
    tested_MLE <- DAISIE:::DAISIE_ML2(
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
      tol = c(0.01, 0.1, 0.001),
      res = 15,
      tolint = c(0.1, 0.01)
    )
  ))
  expected_MLE <- data.frame(
    lambda_c = c(0.0,
                 0.333180238362478,
                 0.0,
                 0.0),
    mu = c(2.718309102127384,
           2.718309102127384,
           2.718309102127384,
           2.718309102127384),
    K = c(Inf,
          Inf,
          Inf,
          Inf),
    gamma = c(0.1389101810794861,
              0.3240600655694914,
              0.1389101810794861,
              0.1389101810794861),
    lambda_a = c(1.17278055953365,
                 1.17278055953365,
                 1.17278055953365 ,
                 1.17278055953365),
    loglik = c(-411.9377194984208,
               -411.9377194984208,
               -411.9377194984208,
               -411.9377194984208),
    df = c(5L, 5L, 5L, 5L),
    conv = c(0L, 0L, 0L, 0L)
  )
  expect_equal(tested_MLE, expected_MLE)
})


test_that("abuse", {
  skip_if(Sys.getenv("CI") == "", message = "Run only on CI")
  expect_error(tested_MLE <- DAISIE:::DAISIE_ML2(
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
