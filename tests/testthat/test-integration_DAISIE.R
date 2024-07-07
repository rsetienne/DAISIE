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
  testthat::expect_equal(loglik, -449.921430187808)
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
  testthat::expect_equal(SR_loglik_CS, loglik_CS)
})

test_that("IW and CS loglik is same when K = Inf", {
  skip_if(Sys.getenv("CI") == "" && !(Sys.getenv("USERNAME") == "rampa"),
          message = "Run only on CI")
  skip_on_cran()
  data(Galapagos_datalist, package = "DAISIE")
  pars1 <- c(0.35, 0.3, Inf, 0.001, 0.3)
  pars2 <- c(80, 11, 1, 0)
  Galapagos_datalist_IW <- Galapagos_datalist
  for(i in 2:9) {
    Galapagos_datalist_IW[[i]]$branching_times <- c(4, 4 - 2*i*0.1,4 -2*i*0.1-0.1)
    Galapagos_datalist_IW[[i]]$stac <- 2
  }

  Galapagos_datalist_IW <- add_brt_table(Galapagos_datalist_IW)
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

  testthat::expect_equal(loglik_IW, loglik_IW2, tol = 5E-6)
  testthat::expect_equal(loglik_IW, loglik_CS, tol = 5E-6)
})

test_that("IW loglik is correct", {
  skip_if(Sys.getenv("CI") == "" && !(Sys.getenv("USERNAME") == "rampa"),
          message = "Run only on CI")
  skip_on_cran()
  frogs_datalist <- NULL
  rm(frogs_datalist)
  data(frogs_datalist, package = "DAISIE")
  frogs_datalist[[1]]$not_present <- frogs_datalist[[1]]$not_present + (M - 300)
  M <- 1000
  ddmodel <- 11
  initparsopt <- c(4.012298e-01,1.699521e-01,1.319595e+02,3.487955e-04)
  idparsopt <- c(1,2,3,4)
  parsfix <- 0
  idparsfix <- 5
  verbose <- 1
  cond <- 1
  res <- 200
  methode <- 'odeint::runge_kutta_fehlberg78'
  tolint <- c(1E-16, 1E-14)

  trparsopt <- initparsopt / (1 + initparsopt)
  trparsopt[which(initparsopt == Inf)] <- 1
  trparsfix <- parsfix / (1 + parsfix)
  trparsfix[which(parsfix == Inf)] <- 1
  pars2 <- c(res, ddmodel, cond, verbose)
  initloglik_IW <- DAISIE_loglik_IW_choosepar(trparsopt = trparsopt,
                                              trparsfix = trparsfix,
                                              idparsopt = idparsopt,
                                              idparsfix = idparsfix,
                                              M = M,
                                              pars2 = pars2,
                                              datalist = frogs_datalist,
                                              methode = methode,
                                              abstolint = tolint[1],
                                              reltolint = tolint[2])

  testthat::expect_equal(initloglik_IW, -215.097677998973, tol = 1E-6)
})

test_that("IW loglik does not error when there are recolonizations", {
  skip_if(Sys.getenv("CI") == "" && !(Sys.getenv("USERNAME") == "rampa"),
          message = "Run only on CI")
  skip_on_cran()
  frogs_sim_CS_T30_M1000_R5000 <- NULL
  rm(frogs_sim_CS_T30_M1000_R5000)
  data(frogs_sim_CS_T30_M1000_R5000, package = "DAISIE")
  datalist <- frogs_sims[[897]]
  M <- 1000
  ddmodel <- 11
  initparsopt <- c(4.012298e-01,1.699521e-01,1.319595e+02,3.487955e-04)
  idparsopt <- c(1,2,3,4)
  parsfix <- 0
  idparsfix <- 5
  verbose <- 1
  cond <- 1
  res <- 200
  methode <- 'odeint::runge_kutta_fehlberg78'
  tolint <- c(1E-16, 1E-14)

  trparsopt <- initparsopt / (1 + initparsopt)
  trparsopt[which(initparsopt == Inf)] <- 1
  trparsfix <- parsfix / (1 + parsfix)
  trparsfix[which(parsfix == Inf)] <- 1
  pars2 <- c(res, ddmodel, cond, verbose)
  initloglik_IW <- DAISIE_loglik_IW_choosepar(trparsopt = trparsopt,
                                              trparsfix = trparsfix,
                                              idparsopt = idparsopt,
                                              idparsfix = idparsfix,
                                              M = M,
                                              pars2 = pars2,
                                              datalist = datalist,
                                              methode = methode,
                                              abstolint = tolint[1],
                                              reltolint = tolint[2])

  testthat::expect_equal(initloglik_IW, -244.7340301475871911, tol = 1E-6)
})

test_that("DAISIE_ML simple case works", {
  skip_if(Sys.getenv("CI") == "" && !(Sys.getenv("USERNAME") == "rampa"),
          message = "Run only on CI")
  skip_on_cran()
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
  testthat::expect_equal(expected_mle, tested_mle)
})

test_that("DAISIE_ML simple case works with zero probability of initial presence", {
  skip_if(Sys.getenv("CI") == "" && !(Sys.getenv("USERNAME") == "rampa"),
          message = "Run only on CI")
  skip_on_cran()
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

    tested_mle <- DAISIE_ML(
      datalist = Galapagos_datalist,
      initparsopt = c(2.5, 2.7, 20, 0.009, 1.01),
      ddmodel = 11,
      idparsopt = 1:5,
      parsfix = 0,
      idparsfix = 6,
      verbose = 0
    )
    expected_calculated_mle <- DAISIE_ML(
      datalist = Galapagos_datalist,
      initparsopt = c(2.5, 2.7, 20, 0.009, 1.01),
      ddmodel = 11,
      idparsopt = 1:5,
      parsfix = NULL,
      idparsfix = NULL,
      verbose = 0
    )

  testthat::expect_equal(expected_mle, tested_mle)
  # Results match if prob_init_pres is removed
  testthat::expect_equal(expected_calculated_mle, tested_mle[-6])
})

test_that("DAISIE_ML simple case works with nonzero probability of initial
          presence", {
            skip_if(Sys.getenv("CI") == "" && !(Sys.getenv("USERNAME") == "rampa"),
                    message = "Run only on CI")
            expected_mle <- data.frame(
              lambda_c = 3.30567366427796,
              mu = 3.86584745010284,
              K = Inf,
              gamma = 0.0144177568387567,
              lambda_a = 0.699608034134341,
              prob_init_pres = 0.1,
              loglik = -78.9245109502749,
              df = 5L,
              conv = 0L
            )
            utils::data(Galapagos_datalist)

            tested_mle <- DAISIE_ML(
              datalist = Galapagos_datalist,
              initparsopt = c(2.5, 2.7, 20, 0.009, 1.01),
              ddmodel = 11,
              idparsopt = 1:5,
              parsfix = 0.1,
              idparsfix = 6,
              verbose = 0
              # verbose = 3,
            )
            testthat::expect_equal(expected_mle, tested_mle, tolerance = 2E-3)
            # tolerance due to different OS results between windows, macOS and
            # ubuntu added in #162
          })


test_that("DAISIE_ML with nonzero probability of initial presence gives
          different results from DAISIE_ML with 0 probability of initial
          presence", {
  skip_if(Sys.getenv("CI") == "" && !(Sys.getenv("USERNAME") == "rampa"),
          message = "Run only on CI")
  utils::data(Galapagos_datalist)

  tested_mle_zero <- DAISIE_ML(
    datalist = Galapagos_datalist,
    initparsopt = c(2.5, 2.7, 20, 0.009, 1.01),
    ddmodel = 11,
    idparsopt = 1:5,
    parsfix = 0,
    idparsfix = 6,
    verbose = 0
  )
  tested_mle_nonzero <- DAISIE_ML(
    datalist = Galapagos_datalist,
    initparsopt = c(2.5, 2.7, 20, 0.009, 1.01),
    ddmodel = 11,
    idparsopt = 1:5,
    parsfix = 0.1,
    idparsfix = 6,
    verbose = 0
  )
  testthat::expect_false(isTRUE(all.equal(tested_mle_zero, tested_mle_nonzero)))
})

test_that("DAISIE_ML gives a -Inf loglikelhood when probability of initial
          presence is nonzero and colonization rate equals 0 for Galapagos", {
            skip_if(Sys.getenv("CI") == "" && !(Sys.getenv("USERNAME") == "rampa"),
                    message = "Run only on CI")
            skip_on_cran()

            utils::data(Galapagos_datalist)
            tested_mle <- DAISIE_ML(
                datalist = Galapagos_datalist,
                initparsopt = c(2.5, 2.7, 20, 0, 1.01, 0.001),
                ddmodel = 11,
                idparsopt = 1:6,
                parsfix = NULL,
                idparsfix = NULL
            )
            testthat::expect_true(is.na(tested_mle$loglik))
          })


test_that("DAISIE_ML simple case works with estimating probability of initial
          presence", {
            skip_if(Sys.getenv("CI") == "" && !(Sys.getenv("USERNAME") == "rampa"),
                    message = "Run only on CI")
            skip_on_cran()

            expected_mle <- data.frame(
              lambda_c = 2.54079308283855,
              mu = 2.66563367593515,
              K = 6249.71023359369,
              gamma = 0.00919247416324124,
              lambda_a = 1.01076206116211,
              prob_init_pres = 9.45796543536632e-06,
              loglik = -75.9935681347126,
              df = 6L,
              conv = 0L
            )

            utils::data(Galapagos_datalist)
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
            testthat::expect_equal(tested_mle, expected_mle)
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
  testthat::expect_equal(fit$loglik, -74.7557, tol = 1E-3)
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
  testthat::expect_equal(loglik_CS_2type_cond0, -61.7028188767349)

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
  testthat::expect_equal(loglik_CS_1type_cond1, -96.45757823017264)

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
  testthat::expect_equal(loglik_CS_2type_cond1, -61.4375954386635)

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
  testthat::expect_equal(loglik_CS_1type_cond5, -95.14000237210625)

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
  testthat::expect_equal(loglik_CS_2type_cond5, -61.3735194058527)
})

test_that("ML with fixed parameters should be different from free parameters
          and be nonzero", {
  skip_if(Sys.getenv("CI") == "" && !(Sys.getenv("USERNAME") == "rampa"),
          message = "Run only on CI")
  skip_on_cran()
  utils::data(Galapagos_datalist)

  tested_mle_free <- DAISIE_ML(
    datalist = Galapagos_datalist,
    initparsopt = c(2.5, 2.7, 20, 0.009, 1.01),
    ddmodel = 11,
    idparsopt = 1:5,
    parsfix = NULL,
    idparsfix = NULL,
    tol = c(1e-2, 1e-3, 1e-4),
    verbose = 0
  )
  tested_mle_fix_clado <- DAISIE_ML(
    datalist = Galapagos_datalist,
    initparsopt = c(2.7, 20, 0.009, 1.01),
    ddmodel = 11,
    idparsopt = 2:5,
    parsfix = 2.5,
    idparsfix = 1,
    tol = c(1e-2, 1e-3, 1e-4),
    verbose = 0
  )
  tested_mle_fix_mu <- DAISIE_ML(
    datalist = Galapagos_datalist,
    initparsopt = c(2.5, 20, 0.009, 1.01),
    ddmodel = 11,
    idparsopt = c(1, 3:5),
    parsfix = 2.7,
    idparsfix = 2,
    tol = c(1e-2, 1e-3, 1e-4),
    verbose = 0
  )
  tested_mle_fix_k <- DAISIE_ML(
    datalist = Galapagos_datalist,
    initparsopt = c(2.5, 2.7, 0.009, 1.01),
    ddmodel = 11,
    idparsopt = c(1, 2, 4, 5),
    parsfix = 20,
    idparsfix = 3,
    tol = c(1e-2, 1e-3, 1e-4),
    verbose = 0
  )
  tested_mle_fix_immig <- DAISIE_ML(
    datalist = Galapagos_datalist,
    initparsopt = c(2.5, 2.7, 20, 1.01),
    ddmodel = 11,
    idparsopt = c(1:3, 5),
    parsfix = 0.009,
    idparsfix = 4,
    tol = c(1e-2, 1e-3, 1e-4),
    verbose = 0
  )
  tested_mle_fix_ana <- DAISIE_ML(
    datalist = Galapagos_datalist,
    initparsopt = c(2.5, 2.7, 20, 0.009),
    ddmodel = 11,
    idparsopt = 1:4,
    parsfix = 1.01,
    idparsfix = 5,
    tol = c(1e-2, 1e-3, 1e-4),
    verbose = 0
  )

  # Fixing one parameter should not return the same as a leaving all free
  testthat::expect_false(isTRUE(all.equal(tested_mle_free, tested_mle_fix_clado)))
  testthat::expect_false(isTRUE(all.equal(tested_mle_free, tested_mle_fix_mu)))
  testthat::expect_false(isTRUE(all.equal(tested_mle_free, tested_mle_fix_k)))
  testthat::expect_false(isTRUE(all.equal(tested_mle_free, tested_mle_fix_immig)))
  testthat::expect_false(isTRUE(all.equal(tested_mle_free, tested_mle_fix_ana)))


  # Fixing one parameter should not return pars as 0, unless set to it
  testthat::expect_false(0 %in% tested_mle_fix_clado[1:5])
  testthat::expect_false(0 %in% tested_mle_fix_mu[1:5])
  testthat::expect_false(0 %in% tested_mle_fix_k[1:5])
  testthat::expect_false(0 %in% tested_mle_fix_immig[1:5])
  testthat::expect_false(0 %in% tested_mle_fix_ana[1:5])
})
