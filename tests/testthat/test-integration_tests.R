context("integration test")
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
  if (Sys.getenv("TRAVIS") != "" || Sys.getenv("APPVEYOR") != "") {
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
  } else {
    testthat::skip("Run only on Travis or AppVeyor")
  }
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
  if (Sys.getenv("TRAVIS") != "" || Sys.getenv("APPVEYOR") != "") {
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
    utils::data(Galapagos_datalist, package = "DAISIE")
    tested_mle <- DAISIE::DAISIE_ML(
      datalist = Galapagos_datalist,
      initparsopt = c(2.5, 2.7, 20, 0.009, 1.01),
      ddmodel = 11,
      idparsopt = 1:5,
      parsfix = NULL,
      idparsfix = NULL
    )
    testthat::expect_equal(expected_mle, tested_mle)
  } else {
    testthat::skip("Run only on Travis or AppVeyor")
  }
})

test_that("The parameter choice for 2type DAISIE_ML works", {
  Galapagos_datalist_2types <- NULL
  rm(Galapagos_datalist_2types)
  utils::data(Galapagos_datalist_2types, package = "DAISIE")
  set.seed(1)
  # MLE and high tolerance for speed-up
  fit <- DAISIE::DAISIE_ML(
    datalist = Galapagos_datalist_2types,
    initparsopt = c(2.183336, 2.517413, 0.009909, 1.080458, 1.316296, 0.001416),
    idparsopt = c(1, 2, 4, 5, 7, 11),
    parsfix = c(Inf, Inf),
    idparsfix = c(3, 8),
    idparsnoshift = c(6, 9, 10),
    res = 30,
    tol = c(1, 1, 1),
    maxiter = 30,
    verbose = FALSE
  )
  testthat::expect_equal(fit$loglik, -74.7557, tol = 1E-3)
})

test_that("DAISIE_sim ontogeny integration", {
  if (Sys.getenv("TRAVIS") != "" || Sys.getenv("APPVEYOR") != "") {
    skip("Temporary skip")
    set.seed(1)
    n_mainland_species <- 50
    island_age <- 5
    clado_rate <- 0.001 # cladogenesis rate
    ext_rate <- 2.683454548 # extinction rate (not used)
    clade_carr_cap <- 0.05  # clade-level carrying capacity
    imm_rate <- 0.01 # immigration rate
    ana_rate <- 0.1 # anagenesis rate
    replicates <- 1
    max_area <- 1000
    peak_time <- 0.1
    sharpness <- 1
    total_island_age <- 15
    sea_level_amplitude <- 0
    sea_level_frequency <- 0
    island_gradient_angle <- 0
    mu_min <- 0.05
    mu_max <- 20
    island_ontogeny <- "beta"
    sea_level <- "const"
    area_pars <- DAISIE::create_area_pars(
      max_area,
      peak_time,
      sharpness,
      total_island_age = total_island_age,
      sea_level_amplitude = sea_level_amplitude,
      sea_level_frequency = sea_level_frequency,
      island_gradient_angle = island_gradient_angle
    )

    expect_silent(
      out <- DAISIE::DAISIE_sim_time_dependent(
        time = island_age,
        M = n_mainland_species,
        pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
        replicates = replicates,
        island_ontogeny = island_ontogeny,
        sea_level = sea_level,
        area_pars = area_pars,
        ext_pars = c(mu_min, mu_max),
        plot_sims = FALSE,
        verbose = FALSE
      )
    )

    expected <- data.frame(
      lambda_c = 0.50851079075334016,
      mu = 2.4663928323693416e-06,
      K = 2.9393253437755784,
      gamma = 0.0079838040792006883,
      lambda_a = 1.3985419762245664e-05,
      loglik = -15.009196369980977,
      df = 5L,
      conv = 0L
    )
    tested <- DAISIE::DAISIE_ML_CS(
      datalist = out[[1]],
      initparsopt = c(2.183336, 2.517413, 20, 1.080458, 1.316296),
      idparsopt = 1:5,
      parsfix = NULL,
      idparsfix = NULL,
      verbose = FALSE
    )
    testthat::expect_equal(expected, tested)
  } else {
    testthat::skip("Run only on Travis or AppVeyor")
  }
})
