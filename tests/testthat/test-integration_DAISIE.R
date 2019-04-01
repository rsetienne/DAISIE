context("integration test")
test_that("loglik simple case works", {
  Galapagos_datalist = NULL
  rm(Galapagos_datalist)
  Galapagos_datalist_2types = NULL
  rm(Galapagos_datalist_2types)
  Macaronesia_datalist = NULL
  rm(Macaronesia_datalist)
  utils::data(Galapagos_datalist_2types, package = "DAISIE")
  pars1 = c(
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
  pars2 = c(100, 11, 0, 0)
  loglik = DAISIE_loglik_all(pars1, pars2, Galapagos_datalist_2types)
  testthat::expect_equal(loglik, -61.7094829913735978)
})

test_that("loglik macaronesia 2 type works", {
  utils::data(Macaronesia_datalist, package = "DAISIE")
  background = c(0, 1.053151832, Inf, 0.052148979, 0.512939011)
  Canaries = c(0.133766934, 1.053151832, Inf, 0.152763179, 0.512939011)
  pars1 = rbind(background, Canaries, background, background)
  pars2 = c(100, 0, 0, 0)
  loglik = 0
  for (i in 1:length(Macaronesia_datalist))
  {
    loglik = loglik + DAISIE_loglik_all(pars1[i, ], pars2, Macaronesia_datalist[[i]], methode = "lsodes")
  }
  testthat::expect_equal(loglik, -454.9347833283220552)
})

test_that("clade specific rate-shift loglik works", {
  utils::data(Galapagos_datalist, package = "DAISIE")
  pars1 = c(0.2, 0.1, Inf, 0.001, 0.3, 0.2, 0.1, Inf, 0.001, 0.3, 1)
  pars2 = c(40, 11, 0, 0)
  SR_loglik_CS = DAISIE_SR_loglik_CS(
    pars1 = pars1,
    pars2 = pars2,
    datalist = Galapagos_datalist,
    methode = 'ode45',
    CS_version = 1
  )
  pars1 = c(0.2, 0.1, Inf, 0.001, 0.3)
  loglik_CS = DAISIE_loglik_CS(
    pars1 = pars1,
    pars2 = pars2,
    datalist = Galapagos_datalist,
    methode = 'ode45',
    CS_version = 1
  )
  testthat::expect_equal(SR_loglik_CS, loglik_CS)
})

test_that("IW and CS loglik is same when K = Inf", {
  utils::data(Galapagos_datalist, package = "DAISIE")
  pars1 = c(0.2, 0.1, Inf, 0.001, 0.3)
  pars2 = c(40, 11, 0, 0)
  loglik_IW = DAISIE_loglik_IW(
    pars1 = pars1,
    pars2 = pars2,
    datalist = Galapagos_datalist,
    methode = 'ode45'
  )
  loglik_CS = DAISIE_loglik_CS(
    pars1 = pars1,
    pars2 = pars2,
    datalist = Galapagos_datalist,
    methode = 'ode45',
    CS_version = 1
  )
  testthat::expect_lt(abs(loglik_IW - loglik_CS), 5E-6)
})

test_that("ontogeny and null-ontogeny loglik is same
          when ontogeny is constant", {
            pars1 = c(0.2, 0.1, 17, 0.001, 0.3)
            pars2 = c(40, 11, 0, 0)
            loglik_CS <- DAISIE_loglik_all(
              pars1 = pars1,
              pars2 = pars2,
              datalist = Galapagos_datalist,
              methode = 'ode45'
            )
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
            pars2 <- c(pars2, translate_island_ontogeny('const'))
            loglik_time <- DAISIE_loglik_all(
              pars1 = pars1_td,
              pars2 = pars2,
              datalist = Galapagos_datalist,
              methode = "ode45"
            )
            testthat::expect_equal(loglik_time, loglik_CS)
})

testthat::test_that("DAISIE_ML simple case works", {
  
  
  if (Sys.getenv("TRAVIS") != "") {
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
  tested_mle <- DAISIE_ML(
    datalist = Galapagos_datalist,
    initparsopt = c(2.5,2.7,20,0.009,1.01),
    ddmodel = 11,
    idparsopt = 1:5,
    parsfix = NULL,
    idparsfix = NULL
  )
  testthat::expect_equal(expected_mle, tested_mle)
  } else {
    skip("Run only on Travis")
  }
})