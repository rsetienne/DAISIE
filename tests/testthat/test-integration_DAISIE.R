test_that("loglik Galapagos works", {
  Galapagos_datalist <- NULL
  rm(Galapagos_datalist)
  Galapagos_datalist_2types <- NULL
  rm(Galapagos_datalist_2types)
  utils::data(Galapagos_datalist_2types, package = "DAISIE")
  pars1 <- c(
    0.195442017,
    0.087959583,
    Inf,
    0.002247364,
    0.873605049,
    3755.202241,
    8.909285094,
    14.99999923,
    0.002247364,
    0.873605049,
    0.163
  )
  pars2 <- c(100, 11, 0, 0)
  loglik <- DAISIE::DAISIE_loglik_all(pars1, pars2, Galapagos_datalist_2types)
  testthat::expect_equal(loglik, -61.7094829913735978)
})

test_that("loglik macaronesia 2 type works", {
  Macaronesia_datalist <- NULL
  rm(Macaronesia_datalist)
  utils::data(Macaronesia_datalist, package = "DAISIE")
  background <- c(0, 1.053151832, Inf, 0.052148979, 0.512939011)
  Canaries <- c(0.133766934, 1.053151832, Inf, 0.152763179, 0.512939011)
  pars1 <- rbind(background, Canaries, background, background)
  pars2 <- c(100, 0, 0, 0)
  loglik <- 0
  for (i in seq_along(Macaronesia_datalist)) {
    loglik <- loglik + DAISIE::DAISIE_loglik_all(pars1[i, ],
                                                 pars2,
                                                 Macaronesia_datalist[[i]],
                                                 methode = "lsodes")
  }
  testthat::expect_equal(loglik, -454.9347833283220552)
})

test_that("clade specific rate-shift loglik works", {
  utils::data(Galapagos_datalist, package = "DAISIE")
  pars1 <- c(0.2, 0.1, Inf, 0.001, 0.3, 0.2, 0.1, Inf, 0.001, 0.3, 1)
  pars2 <- c(40, 11, 0, 0)
  SR_loglik_CS <- DAISIE::DAISIE_SR_loglik_CS(
    pars1 = pars1,
    pars2 = pars2,
    datalist = Galapagos_datalist,
    methode = "ode45",
    CS_version = 1)
  pars1 <- c(0.2, 0.1, Inf, 0.001, 0.3)
  loglik_CS <- DAISIE::DAISIE_loglik_CS(
    pars1 = pars1,
    pars2 = pars2,
    datalist = Galapagos_datalist,
    methode = "ode45",
    CS_version = 1)
  testthat::expect_equal(SR_loglik_CS, loglik_CS)
})

test_that("IW and CS loglik is same when K = Inf", {
  skip_if(Sys.getenv("CI") == "" || !(Sys.getenv("USERNAME") == "rampa"),
          message = "Run only on CI")
  utils::data(Galapagos_datalist, package = "DAISIE")
  pars1 <- c(0.2, 0.1, Inf, 0.001, 0.3)
  pars2 <- c(40, 11, 0, 0)
  loglik_IW <- DAISIE::DAISIE_loglik_IW(
    pars1 = pars1,
    pars2 = pars2,
    datalist = Galapagos_datalist,
    methode = "ode45")
  loglik_CS <- DAISIE::DAISIE_loglik_CS(
    pars1 = pars1,
    pars2 = pars2,
    datalist = Galapagos_datalist,
    methode = "ode45",
    CS_version = 1)
  testthat::expect_lt(abs(loglik_IW - loglik_CS), 5E-6)
})

test_that("ontogeny and null-ontogeny loglik is same when ontogeny is
          constant", {
            skip("Temporary skip")
            pars1 <- c(0.2, 0.1, 17, 0.001, 0.3)
            pars2 <- c(40, 11, 0, 0)
            utils::data(Galapagos_datalist, package = "DAISIE")
            loglik_CS <- DAISIE::DAISIE_loglik_all(
              pars1 = pars1,
              pars2 = pars2,
              datalist = Galapagos_datalist,
              methode = "ode45")
            pars1_td <- c(
              max_area = 1,
              proportional_peak_t = 0.2,
              peak_sharpness = 1,
              total_island_age = 15,
              lac = pars1[1],
              mu_min = pars1[2],
              mu_max = pars1[2],
              K0 = pars1[3],
              gam = pars1[4],
              laa = pars1[5]
            )
            pars1_td <- DAISIE:::order_pars1(pars1_td)
            pars2 <- c(pars2, DAISIE::translate_island_ontogeny("const"))
            loglik_time <- DAISIE::DAISIE_loglik_all(
              pars1 = pars1_td,
              pars2 = pars2,
              datalist = Galapagos_datalist,
              methode = "ode45"
            )
            testthat::expect_equal(loglik_time, loglik_CS)
})

testthat::test_that("DAISIE_ML simple case works", {
  skip_if(Sys.getenv("CI") == "" || !(Sys.getenv("USERNAME") == "rampa"),
          message = "Run only on CI")
  expected_mle <- data.frame(
    lambda_c = 2.55847849219339,
    mu = 2.68768191590176,
    K = 6765.0637400135,
    gamma = 0.00932987953669849,
    lambda_a = 1.00838182578826,
    loglik = -76.0001379108545,
    df = 5L,
    conv = 0L
  )
  utils::data(Galapagos_datalist)
  cat("\n")
  tested_mle <- DAISIE::DAISIE_ML(
    datalist = Galapagos_datalist,
    initparsopt = c(2.5, 2.7, 20, 0.009, 1.01),
    ddmodel = 11,
    idparsopt = 1:5,
    parsfix = NULL,
    idparsfix = NULL
  )
  testthat::expect_equal(expected_mle, tested_mle)
})

test_that("The parameter choice for 2type DAISIE_ML works", {
  Galapagos_datalist_2types <- NULL
  rm(Galapagos_datalist_2types)
  utils::data(Galapagos_datalist_2types, package = "DAISIE")
  set.seed(1)
  # MLE and high tolerance for speed-up
  cat("\n")
  fit <- DAISIE::DAISIE_ML(
    datalist = Galapagos_datalist_2types,
    initparsopt = c(2.183336, 2.517413, 0.009909, 1.080458, 1.316296, 0.001416),
    idparsopt = c(1, 2, 4, 5, 7, 11),
    parsfix = c(Inf, Inf),
    idparsfix = c(3, 8),
    idparsnoshift = c(6, 9, 10),
    res = 30,
    tol = c(1, 1, 1),
    maxiter = 30
  )
  testthat::expect_equal(fit$loglik, -74.7557, tol = 1E-3)
})

test_that("conditioning works", {
  # Cond 0
  ## 1 type
  utils::data(Galapagos_datalist, package = "DAISIE")
  pars1_1type_cond0 <- c(0.2, 0.1, Inf, 0.001, 0.3)
  pars2_1type_cond0 <- c(40, 11, 0, 0)
  loglik_CS_1type_cond0 <- DAISIE::DAISIE_loglik_CS(
    pars1 = pars1_1type_cond0,
    pars2 = pars2_1type_cond0,
    datalist = Galapagos_datalist,
    methode = "ode45",
    CS_version = 1
  )
  testthat::expect_equal(loglik_CS_1type_cond0, -96.49629968062564)

  ## 2 type
  utils::data(Galapagos_datalist_2types, package = "DAISIE")
  pars1_2type_cond0 <- c(
    0.195442017,
    0.087959583,
    Inf,
    0.002247364,
    0.873605049,
    3755.202241,
    8.909285094,
    14.99999923,
    0.002247364,
    0.873605049,
    0.163
  )
  pars2_2type_cond0 <- c(100, 11, 0, 0)
  loglik_CS_2type_cond0 <- DAISIE::DAISIE_loglik_CS(
    pars1_2type_cond0,
    pars2_2type_cond0,
    Galapagos_datalist_2types
  )
  testthat::expect_equal(loglik_CS_2type_cond0, -61.709482984890265)

  # Cond 1
  ## 1 type
  utils::data(Galapagos_datalist, package = "DAISIE")
  pars1_1type_cond1 <- c(0.2, 0.1, Inf, 0.001, 0.3)
  pars2_1type_cond1 <- c(40, 11, 1, 0)
  loglik_CS_1type_cond1 <- DAISIE::DAISIE_loglik_CS(
    pars1 = pars1_1type_cond1,
    pars2 = pars2_1type_cond1,
    datalist = Galapagos_datalist,
    methode = 'ode45',
    CS_version = 1
  )
  testthat::expect_equal(loglik_CS_1type_cond1, -96.463184608046333)

  ## 2 type
  utils::data(Galapagos_datalist_2types, package = "DAISIE")
  pars1_2type_cond1 <- c(
    0.195442017,
    0.087959583,
    Inf,
    0.002247364,
    0.873605049,
    3755.202241,
    8.909285094,
    14.99999923,
    0.002247364,
    0.873605049,
    0.163
  )
  pars2_2type_cond1 <- c(100, 11, 1, 0)
  loglik_CS_2type_cond1 <- DAISIE::DAISIE_loglik_CS(
    pars1_2type_cond1,
    pars2_2type_cond1,
    Galapagos_datalist_2types
  )
  testthat::expect_equal(loglik_CS_2type_cond1,-61.709153802942346)

  # Cond 5
  ## 1 type
  utils::data(Galapagos_datalist, package = "DAISIE")
  pars1_1type_cond5 <- c(0.2, 0.1, Inf, 0.001, 0.3)
  pars2_1type_cond5 <- c(40, 11, 5, 0)
  loglik_CS_1type_cond5 <- DAISIE::DAISIE_loglik_CS(
    pars1 = pars1_1type_cond5,
    pars2 = pars2_1type_cond5,
    datalist = Galapagos_datalist,
    methode = 'ode45',
    CS_version = 1
  )
  testthat::expect_equal(loglik_CS_1type_cond5,-95.1456087499799423)

  ## 2 type
  utils::data(Galapagos_datalist_2types, package = "DAISIE")
  pars1_2type_cond5 <- c(
    0.195442017,
    0.087959583,
    Inf,
    0.002247364,
    0.873605049,
    3755.202241,
    8.909285094,
    14.99999923,
    0.002247364,
    0.873605049,
    0.163
  )
  pars2_2type_cond5 <- c(100, 11, 5, 0)
  loglik_CS_2type_cond5 <- DAISIE::DAISIE_loglik_all(
    pars1_2type_cond5,
    pars2_2type_cond5,
    Galapagos_datalist_2types
  )
  testthat::expect_equal(loglik_CS_2type_cond5, -61.5667762281177673)
})

test_that("various solver options give similar results", {
  # Test is not included in coverage due to issue with running loglik_IW
  # code from covr::package_coverage()

  utils::data(frogs_datalist, package = "DAISIE") #nocov start
  pars1 <- c(0.2, 0.1, 1000.1, 0.001, 0.3)
  pars2 <- c(40, 11, 0, 0)

  IW0 <- DAISIE::DAISIE_loglik_IW(
    pars1 = pars1,
    pars2 = pars2,
    datalist = frogs_datalist,
    methode = 'odeint::runge_kutta_fehlberg78',
    abstolint = 1E-12,
    reltolint = 1E-10,
  )
  IW1 <- DAISIE::DAISIE_loglik_IW(
    pars1 = pars1,
    pars2 = pars2,
    datalist = frogs_datalist,
    methode = 'odeint::runge_kutta_cash_karp54',
    abstolint = 1E-12,
    reltolint = 1E-10,
  )
  IW2 <- DAISIE::DAISIE_loglik_IW(
    pars1 = pars1,
    pars2 = pars2,
    datalist = frogs_datalist,
    methode = 'odeint::runge_kutta_dopri5',
    abstolint = 1E-12,
    reltolint = 1E-10,
  )
  IW3 <- DAISIE::DAISIE_loglik_IW(
    pars1 = pars1,
    pars2 = pars2,
    datalist = frogs_datalist,
    methode = 'odeint::bulirsch_stoer',
    abstolint = 1E-12,
    reltolint = 1E-10,
  )
  IW4 <- DAISIE::DAISIE_loglik_IW(
    pars1 = pars1,
    pars2 = pars2,
    datalist = frogs_datalist,
    methode = 'ode45',
    abstolint = 1E-12,
    reltolint = 1E-10,
  )

  #print(c(IW0,IW1,IW2,IW3,IW4))
  testthat::expect_equal(IW0,IW1, tolerance = 1E-4)
  testthat::expect_equal(IW0,IW2, tolerance = 1E-4)
  testthat::expect_equal(IW0,IW3, tolerance = 1E-4)
  testthat::expect_equal(IW0,IW4, tolerance = 1E-4) #nocov end
})
