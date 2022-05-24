test_that("loglik Galapagos works", {
  Galapagos_datalist <- NULL
  rm(Galapagos_datalist)
  Galapagos_datalist_2types <- NULL
  rm(Galapagos_datalist_2types)
  data(Galapagos_datalist_2types, package = "DAISIE")
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
  loglik <- DAISIE_loglik_all(pars1, pars2, Galapagos_datalist_2types)
  testthat::expect_equal(loglik, -61.70281911731144)
})

test_that("loglik macaronesia 2 type works", {
  Macaronesia_datalist <- NULL
  rm(Macaronesia_datalist)
  data(Macaronesia_datalist, package = "DAISIE")
  background <- c(0, 1.053151832, Inf, 0.052148979, 0.512939011)
  Canaries <- c(0.133766934, 1.053151832, Inf, 0.152763179, 0.512939011)
  pars1 <- rbind(background, Canaries, background, background)
  pars2 <- c(100, 0, 0, 0)
  loglik <- 0
  for (i in seq_along(Macaronesia_datalist)) {
    loglik <- loglik + DAISIE_loglik_all(pars1[i, ],
                                                 pars2,
                                                 Macaronesia_datalist[[i]],
                                                 methode = "lsodes")
  }
  expect_equal(loglik, -449.921430187808)
})

test_that("clade specific rate-shift loglik works", {
  data(Galapagos_datalist, package = "DAISIE")
  pars1 <- c(0.2, 0.1, Inf, 0.001, 0.3, 0.2, 0.1, Inf, 0.001, 0.3, 1)
  pars2 <- c(40, 11, 0, 0)
  SR_loglik_CS <- DAISIE_SR_loglik_CS(
    pars1 = pars1,
    pars2 = pars2,
    datalist = Galapagos_datalist,
    methode = "ode45",
    CS_version = 1)
  pars1 <- c(0.2, 0.1, Inf, 0.001, 0.3)
  loglik_CS <- DAISIE_loglik_CS(
    pars1 = pars1,
    pars2 = pars2,
    datalist = Galapagos_datalist,
    methode = "ode45",
    CS_version = 1)
  expect_equal(SR_loglik_CS, loglik_CS)
})

test_that("IW and CS loglik is same when K = Inf", {
  skip_if(Sys.getenv("CI") == "" && !(Sys.getenv("USERNAME") == "rampa"),
          message = "Run only on CI")
  data(Galapagos_datalist, package = "DAISIE")
  pars1 <- c(0.35, 0.3, Inf, 0.001, 0.3)
  pars2 <- c(80, 11, 1, 0)
  Galapagos_datalist_IW <- Galapagos_datalist
  for(i in 2:9) {
    Galapagos_datalist_IW[[i]]$branching_times <- c(4, 4 - 2*i*0.1,4 -2*i*0.1-0.1)
    Galapagos_datalist_IW[[i]]$stac <- 2
  }

  Galapagos_datalist_IW <- DAISIE:::add_brt_table(Galapagos_datalist_IW)
  loglik_IW <- DAISIE_loglik_IW(
      pars1 = pars1,
      pars2 = pars2,
      datalist = Galapagos_datalist_IW,
      methode = "ode45"
  )

    loglik_IW2 <- DAISIE_loglik_IW(
      pars1 = pars1,
      pars2 = pars2,
      datalist = Galapagos_datalist_IW,
      methode = "odeint::runge_kutta_fehlberg78"
    )


    loglik_CS <- DAISIE_loglik_CS(
      pars1 = pars1,
      pars2 = pars2,
      datalist = Galapagos_datalist_IW,
      methode = "ode45",
      CS_version = 1
    )

  expect_equal(loglik_IW, loglik_IW2, tol = 5E-6)
  expect_equal(loglik_IW, loglik_CS, tol = 5E-6)
})

test_that("DAISIE_ML simple case works", {
  skip_if(Sys.getenv("CI") == "" && !(Sys.getenv("USERNAME") == "rampa"),
          message = "Run only on CI")
  expected_mle <- data.frame(
    lambda_c = 2.583731356303842,
    mu = 2.708828027514834,
    K = 2992.207701921788,
    gamma = 0.00937711049761019,
    lambda_a = 0.9993246958280274,
    loglik = -75.99266304738612,
    df = 5L,
    conv = 0L
  )
  utils::data(Galapagos_datalist)

  invisible(capture.output(
    tested_mle <- DAISIE_ML(
      datalist = Galapagos_datalist,
      initparsopt = c(2.5, 2.7, 20, 0.009, 1.01),
      ddmodel = 11,
      idparsopt = 1:5,
      parsfix = NULL,
      idparsfix = NULL
    )
  ))
  expect_equal(expected_mle, tested_mle)
})

test_that("DAISIE_ML simple case works with zero probability of initial presence", {
  skip_if(Sys.getenv("CI") == "" && !(Sys.getenv("USERNAME") == "rampa"),
          message = "Run only on CI")
  expected_mle <- data.frame(
    lambda_c = 2.583731356303842,
    mu = 2.708828027514834,
    K = 2992.207701921788,
    gamma = 0.00937711049761019,
    lambda_a = 0.9993246958280274,
    prob_init_pres = 0,
    loglik = -75.99266304738612,
    df = 5L,
    conv = 0L
  )
  utils::data(Galapagos_datalist)

  invisible(capture.output(
    tested_mle <- DAISIE_ML(
      datalist = Galapagos_datalist,
      initparsopt = c(2.5, 2.7, 20, 0.009, 1.01),
      ddmodel = 11,
      idparsopt = 1:5,
      parsfix = 0,
      idparsfix = 6
    )
  ))
  expect_equal(expected_mle, tested_mle)
})

test_that("DAISIE_ML simple case works with nonzero probability of initial presence", {
  skip_if(Sys.getenv("CI") == "" && !(Sys.getenv("USERNAME") == "rampa"),
          message = "Run only on CI")
  expected_mle <- data.frame(
    lambda_c = 2.58373135630384,
    mu = 2.70882802751483,
    K = 2992.20770192179,
    gamma = 0.00937711049761019,
    lambda_a = 0.999324695828027,
    prob_init_pres = 0.1,
    loglik = -75.9926628720867,
    df = 5L,
    conv = 0L
  )
  utils::data(Galapagos_datalist)

  invisible(capture.output(
    tested_mle <- DAISIE_ML(
      datalist = Galapagos_datalist,
      initparsopt = c(2.5, 2.7, 20, 0.009, 1.01),
      ddmodel = 11,
      idparsopt = 1:5,
      parsfix = 0.1,
      idparsfix = 6
    )
  ))
  expect_equal(expected_mle, tested_mle)
})

test_that("DAISIE_ML simple case works with estimating probability of initial presence", {
  skip_if(Sys.getenv("CI") == "" && !(Sys.getenv("USERNAME") == "rampa"),
          message = "Run only on CI")

  if (identical(Sys.getenv("OS"), "Windows_NT")) {
    expected_mle <- data.frame(
      lambda_c = 2.53430497145461,
      mu = 2.66658569091753,
      K = 2136343.97554965,
      gamma = 0.00930345848936764,
      lambda_a = 1.0119011474385,
      prob_init_pres = 3.21939792431987e-10,
      loglik = -75.9925548510873,
      df = 6L,
      conv = 0L
    )
  } else {
    expected_mle <- data.frame(
      lambda_c = 2.53429041285525,
      mu = 2.66553367929804,
      K = 3876287.99373951,
      gamma = 0.00929455817164771,
      lambda_a = 1.01208298276806,
      prob_init_pres = 1.39803679789886e-08,
      loglik = -75.992565711427,
      df = 6L,
      conv = 0L
    )
  }


  utils::data(Galapagos_datalist)
options(digits = 15)
  invisible(capture.output(
    tested_mle <- DAISIE_ML(
      datalist = Galapagos_datalist,
      initparsopt = c(2.5, 2.7, 20, 0.009, 1.01, 0.001),
      ddmodel = 11,
      idparsopt = 1:6,
      parsfix = NULL,
      idparsfix = NULL
    )
  ))
  print(tested_mle)
  expect_equal(tested_mle, expected_mle)
})

test_that("The parameter choice for 2type DAISIE_ML works", {
  Galapagos_datalist_2types <- NULL
  rm(Galapagos_datalist_2types)
  data(Galapagos_datalist_2types, package = "DAISIE")
  set.seed(1)
  # MLE and high tolerance for speed-up

  invisible(capture.output(
    fit <- DAISIE_ML(
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
  ))
  expect_equal(fit$loglik, -74.7557, tol = 1E-3)
})

test_that("conditioning works", {
  # Cond 0
  ## 1 type
  data(Galapagos_datalist, package = "DAISIE")
  pars1_1type_cond0 <- c(0.2, 0.1, Inf, 0.001, 0.3)
  #pars1_1type_cond0 <- c(0, 0, Inf, 1, 0)
  pars2_1type_cond0 <- c(40, 11, 0, 0)
  res1 <- DAISIE_loglik_CS(
    pars1 = pars1_1type_cond0,
    pars2 = pars2_1type_cond0,
    datalist = Galapagos_datalist,
    methode = "ode45",
    CS_version = 1
  )
  res2 <- DAISIE_loglik_CS(
    pars1 = pars1_1type_cond0,
    pars2 = pars2_1type_cond0,
    datalist = Galapagos_datalist,
    methode = "deSolve_R::ode45",
    CS_version = 1
  )
  res3 <- loglik_CS_1type_cond0 <- DAISIE_loglik_CS(
    pars1 = pars1_1type_cond0,
    pars2 = pars2_1type_cond0,
    datalist = Galapagos_datalist,
    methode = "odeint::runge_kutta_fehlberg78",
    CS_version = 1
  )

  testthat::expect_equal(res1, res3)
  testthat::expect_equal(res2, res3, tol = 1E-4)
  testthat::expect_equal(loglik_CS_1type_cond0, -96.49069330275196)

#  Status of colonist: 0, Parameters: 0.200000 0.100000 Inf 0.001000 0.300000 , Loglikelihood: -0.003424
#  Status of colonist: 1, Parameters: 0.200000 0.100000 Inf 0.001000 0.300000 , Loglikelihood: -6.494398
#  Status of colonist: 4, Parameters: 0.200000 0.100000 Inf 0.001000 0.300000 , Loglikelihood: -7.113751
#  Status of colonist: 2, Parameters: 0.200000 0.100000 Inf 0.001000 0.300000 , Loglikelihood: -31.251817
#  Status of colonist: 2, Parameters: 0.200000 0.100000 Inf 0.001000 0.300000 , Loglikelihood: -14.421388
#  Status of colonist: 2, Parameters: 0.200000 0.100000 Inf 0.001000 0.300000 , Loglikelihood: -8.594293
#  Status of colonist: 2, Parameters: 0.200000 0.100000 Inf 0.001000 0.300000 , Loglikelihood: -10.599996
#  Status of colonist: 1, Parameters: 0.200000 0.100000 Inf 0.001000 0.300000 , Loglikelihood: -6.494398
#  Status of colonist: 2, Parameters: 0.200000 0.100000 Inf 0.001000 0.300000 , Loglikelihood: -8.123768
#  Parameters: 0.200000 0.100000 Inf 0.001000 0.300000 , Loglikelihood: -96.490693

  ## 2 type
  data(Galapagos_datalist_2types, package = "DAISIE")
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
  loglik_CS_2type_cond0 <- DAISIE_loglik_CS(
    pars1_2type_cond0,
    pars2_2type_cond0,
    Galapagos_datalist_2types
  )
  expect_equal(loglik_CS_2type_cond0, -61.7028188767349)

  # Cond 1
  ## 1 type
  data(Galapagos_datalist, package = "DAISIE")
  pars1_1type_cond1 <- c(0.2, 0.1, Inf, 0.001, 0.3)
  pars2_1type_cond1 <- c(40, 11, 1, 0)
  loglik_CS_1type_cond1 <- DAISIE_loglik_CS(
    pars1 = pars1_1type_cond1,
    pars2 = pars2_1type_cond1,
    datalist = Galapagos_datalist,
    methode = 'ode45',
    CS_version = 1
  )
  expect_equal(loglik_CS_1type_cond1, -96.45757823017264)

  ## 2 type
  data(Galapagos_datalist_2types, package = "DAISIE")
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
  loglik_CS_2type_cond1 <- DAISIE_loglik_CS(
    pars1_2type_cond1,
    pars2_2type_cond1,
    Galapagos_datalist_2types
  )
  expect_equal(loglik_CS_2type_cond1, -61.4375954386635)

  # Cond 5
  ## 1 type
  utils::data(Galapagos_datalist, package = "DAISIE")
  pars1_1type_cond5 <- c(0.2, 0.1, Inf, 0.001, 0.3)
  pars2_1type_cond5 <- c(40, 11, 5, 0)
  loglik_CS_1type_cond5 <- DAISIE_loglik_CS(
    pars1 = pars1_1type_cond5,
    pars2 = pars2_1type_cond5,
    datalist = Galapagos_datalist,
    methode = 'ode45',
    CS_version = 1
  )
  expect_equal(loglik_CS_1type_cond5, -95.14000237210625)

  ## 2 type
  data(Galapagos_datalist_2types, package = "DAISIE")
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
  loglik_CS_2type_cond5 <- DAISIE_loglik_all(
    pars1_2type_cond5,
    pars2_2type_cond5,
    Galapagos_datalist_2types
  )
  expect_equal(loglik_CS_2type_cond5, -61.3735194058527)
})
