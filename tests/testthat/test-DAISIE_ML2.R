test_that("use", {
  skip_if(Sys.getenv("CI") == "", message = "Run only on CI")
  skip_on_cran()
  utils::data(Macaronesia_datalist, package = "DAISIE")
  invisible(capture.output(
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
      tol = c(0.01, 0.1, 0.001),
      res = 15,
      tolint = c(0.1, 0.01)
    )
  ))
  expected_MLE <- data.frame(
    lambda_c = c(0.0,
                 0.315015693879739,
                 0.0,
                 0.0),
    mu = c(2.67928476648249,
           2.67928476648249,
           2.67928476648249,
           2.67928476648249),
    K = c(Inf,
          Inf,
          Inf,
          Inf),
    gamma = c(0.141748213992604,
              0.338044376412124,
              0.141748213992604,
              0.141748213992604),
    lambda_a = c(1.10809210249609,
                 1.10809210249609,
                 1.10809210249609 ,
                 1.10809210249609),
    loglik = c(-409.60117813893,
               -409.60117813893,
               -409.60117813893,
               -409.60117813893),
    df = c(5L, 5L, 5L, 5L),
    conv = c(0L, 0L, 0L, 0L)
  )
  expect_equal(tested_MLE, expected_MLE)
})


test_that("abuse", {
  skip_if(Sys.getenv("CI") == "", message = "Run only on CI")
  skip_on_cran()
  expect_error(tested_MLE <- DAISIE_ML2(
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
