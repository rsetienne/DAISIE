test_that("DAISIE_loglik_CS_choice produces correct output for CS_version 1", {
  #skip_if(Sys.getenv("CI") == "", message = "Run only on CI")
  #skip_on_cran()
  pars1 <- c(2.000, 2.700, 20.000, 0.009, 1.010)
  pars2 <- c(1.0e+02, 1.1e+01, 0.0e+00, 0.0e+00, NA)
  brts <- c(4.0000, 3.0282, 1.3227, 0.8223, 0.4286, 0.3462, 0.2450, 0.0808,
            0.0527, 0.0327, 0.0221, 0.1180, 0.0756, 0.0525, 0.0322, 0.0118)
  stac <- 2
  missnumspec <- 0
  loglik <- DAISIE_loglik_CS_choice(pars1 = pars1,
                                    pars2 = pars2,
                                    brts = brts,
                                    stac = stac,
                                    missnumspec = missnumspec)

  testthat::expect_true(is.numeric(loglik))
  testthat::expect_equal(loglik, -17.6535269566649)

})

test_that("DAISIE_loglik_CS_choice produces correct output for relaxed-rate
          model (CS_version = 2)", {
  #skip_if(Sys.getenv("CI") == "", message = "Run only on CI")
  #skip_on_cran()
  pars1 <- c(2.000, 2.700, 20.000, 0.009, 1.010)
  pars2 <- c(1.0e+02, 1.1e+01, 0.0e+00, 0.0e+00)
  brts <- c(4.0000, 3.0282, 1.3227, 0.8223, 0.4286, 0.3462, 0.2450, 0.0808,
            0.0527, 0.0327, 0.0221, 0.1180, 0.0756, 0.0525, 0.0322, 0.0118)
  stac <- 2
  missnumspec <- 0
  CS_version <- list(model = 2,
                     relaxed_par = "cladogenesis",
                     par_sd = 1,
                     par_upper_bound = Inf)

  invisible(capture.output(loglik <- DAISIE_loglik_CS_choice(pars1 = pars1,
                                                             pars2 = pars2,
                                                             brts = brts,
                                                             stac = stac,
                                                             missnumspec = missnumspec,
                                                             CS_version = CS_version)))
  testthat::expect_true(is.numeric(loglik))
  testthat::expect_equal(loglik, -9.55018422213285)

})

test_that("DAISIE_loglik_CS_choice produces same output for CS_version = 0
          (with M = 1) and CS_version = 1 ", {
  #skip_on_cran()
  #skip_if(Sys.getenv("CI") == "", message = "Run only on CI")
  pars1 <- c(2.000, 2.700, 20.000, 0.009, 1.010)
  pars2 <- c(100, 11, 0, 0)
  brts <- c(4.0000, 3.0282, 1.3227, 0.8223, 0.4286, 0.3462, 0.2450, 0.0808,
            0.0527, 0.0327, 0.0221, 0.1180, 0.0756, 0.0525, 0.0322, 0.0118)
  stac <- 2
  missnumspec <- 0
  CS_version <- 0
  datalist <- list(branching_times = brts, stac = stac)
  loglik0 <- DAISIE_loglik_CS_choice(pars1 = pars1,
                                     pars2 = pars2,
                                     datalist = datalist,
                                     brts = brts,
                                     stac = stac,
                                     missnumspec = missnumspec,
                                     CS_version = CS_version)
  CS_version <- 1
  loglik1 <- DAISIE_loglik_CS_choice(pars1 = pars1,
                                     pars2 = pars2,
                                     brts = brts,
                                     stac = stac,
                                     missnumspec = missnumspec,
                                     CS_version = CS_version)

  testthat::expect_equal(loglik0, loglik1)
})

test_that("DAISIE_loglik_all produces correct output for relaxed-rate model", {
  #skip_if(Sys.getenv("CI") == "", message = "Run only on CI")
  #skip_on_cran()
  utils::data(Galapagos_datalist)

  invisible(capture.output(suppressWarnings(
    loglik1 <- DAISIE_loglik_all(
      pars1 = c(2.55068735, 2.68345455, 10.00000000, 0.00933207, 1.01007312),
      pars2 = c(100, 0, 0, 0, NA),
      datalist = Galapagos_datalist,
      methode = "lsodes",
      CS_version = list(model = 2,
                        relaxed_par = "cladogenesis",
                        par_sd = 1,
                        par_upper_bound = Inf),
      abstolint = 1e-16,
      reltolint = 1e-10
    )
  )))

  invisible(capture.output(suppressWarnings(
    loglik2a <- DAISIE_loglik_all(
      pars1 = c(2.55068735, 2.68345455, 10.00000000, 0.00933207, 1.01007312),
      pars2 = c(100, 0, 0, 0, NA),
      datalist = Galapagos_datalist,
      methode = "lsodes",
      CS_version = list(model = 2,
                        relaxed_par = "cladogenesis",
                        par_sd = 1,
                        par_upper_bound = Inf,
                        integration_method = 'MC',
                        sample_size = 100,
                        seed = 42),
      abstolint = 1e-16,
      reltolint = 1e-10
    )
  )))

  invisible(capture.output(suppressWarnings(
    loglik2b <- DAISIE_loglik_all(
      pars1 = c(2.55068735, 2.68345455, 10.00000000, 0.00933207, 1.01007312),
      pars2 = c(100, 0, 0, 0, NA),
      datalist = Galapagos_datalist,
      methode = "odeint::runge_kutta_cash_karp54",
      CS_version = list(model = 2,
                        relaxed_par = "cladogenesis",
                        par_sd = 1,
                        par_upper_bound = Inf,
                        integration_method = 'MC',
                        sample_size = 100,
                        seed = 42),
      abstolint = 1e-16,
      reltolint = 1e-10
    )
  )))

  invisible(capture.output(suppressWarnings(
    loglik2c <- DAISIE_loglik_all(
      pars1 = c(2.55068735, 2.68345455, 10.00000000, 0.00933207, 1.01007312),
      pars2 = c(100, 0, 0, 0, NA),
      datalist = Galapagos_datalist,
      methode = "odeint::runge_kutta_cash_karp54",
      CS_version = list(model = 2,
                        relaxed_par = "cladogenesis",
                        par_sd = 1,
                        par_upper_bound = Inf,
                        integration_method = 'MC',
                        sample_size = 100,
                        seed = 42,
                        parallel = TRUE,
                        n_cores = 8),
      abstolint = 1e-16,
      reltolint = 1e-10
    )
  )))

  invisible(capture.output(suppressWarnings(
    loglik3 <- DAISIE_loglik_all(
      pars1 = c(2.55068735, 2.68345455, 10.00000000, 0.00933207, 1.01007312),
      pars2 = c(100, 0, 0, 0, NA),
      datalist = Galapagos_datalist,
      methode = "lsodes",
      CS_version = list(model = 2,
                        relaxed_par = "cladogenesis",
                        par_sd = 1,
                        par_upper_bound = Inf,
                        integration_method = 'stratified',
                        sample_size = 100),
      abstolint = 1e-16,
      reltolint = 1e-10
    )
  )))

  testthat::expect_true(is.numeric(loglik1))
  testthat::expect_true(is.numeric(loglik2a))
  testthat::expect_true(is.numeric(loglik2b))
  testthat::expect_true(is.numeric(loglik2c))
  testthat::expect_true(is.numeric(loglik3))
  testthat::expect_equal(loglik1, -77.5030062791617)
  testthat::expect_equal(loglik2a,loglik2b)
  testthat::expect_equal(loglik2b,loglik2c)
  testthat::expect_equal(loglik2b, -77.5, tol = 1E-2)
  testthat::expect_equal(loglik3, -77.5, tol = 1E-2)

})

test_that("DAISIE_loglik produces correct output", {
  #skip_if(Sys.getenv("CI") == "", message = "Run only on CI")
  #skip_on_cran()
  output <- DAISIE_loglik(pars1 = c(2.061154e-09, 2.683455e+00, 1.000000e+01,
                                    9.332070e-03, 1.010073e+00),
                          pars2 = c(100, 0, 0, 0, NA),
                          brts = 4,
                          stac = 0,
                          missnumspec = 0,
                          methode = "lsodes",
                          abstolint = 1E-16,
                          reltolint = 1E-10,
                          verbose = FALSE)
  testthat::expect_equal(output, -0.00347317077256095)
})


test_that("DAISIE_loglik_all produces same output for CS_version 0 and 1 with
          and without conditioning", {
  #skip_if(Sys.getenv("CI") == "", message = "Run only on CI")
  #skip_on_cran()
  utils::data(Galapagos_datalist)
  Galapagos_datalist2 <- Galapagos_datalist
  for(i in 2:9) {
    Galapagos_datalist2[[i]]$branching_times <- c(4, 4 - 2*i*0.1,4 -2*i*0.1-0.1)
    Galapagos_datalist2[[i]]$stac <- 2
  }
  Galapagos_datalist2 <- add_brt_table(Galapagos_datalist2)
  loglik_CS00 <- DAISIE_loglik_all(
    pars1 = c(2.55068735, 2.68345455, 10.00000000, 0.00933207, 1.01007312),
    pars2 = c(100, 11, 0, 0, NA),
    datalist = Galapagos_datalist2,
    methode = "odeint::runge_kutta_fehlberg78",
    CS_version = 0,
    abstolint = 1e-16,
    reltolint = 1e-10)
  loglik_CS10 <- DAISIE_loglik_all(
    pars1 = c(2.55068735, 2.68345455, 10.00000000, 0.00933207, 1.01007312),
    pars2 = c(100, 11, 0, 0, NA),
    datalist = Galapagos_datalist2,
    methode = "ode45",
    CS_version = list(model = 1, function_to_optimize = 'DAISIE'),
    abstolint = 1e-16,
    reltolint = 1e-10)
  testthat::expect_equal(loglik_CS00, loglik_CS10, tol = 5E-6)
  loglik_CS01 <- DAISIE_loglik_all(
    pars1 = c(2.55068735, 2.68345455, 10.00000000, 0.00933207, 1.01007312),
    pars2 = c(100, 11, 1, 0, NA),
    datalist = Galapagos_datalist2,
    methode = "odeint::runge_kutta_fehlberg78",
    CS_version = 0,
    abstolint = 1e-16,
    reltolint = 1e-10)
  loglik_CS11 <- DAISIE_loglik_all(
    pars1 = c(2.55068735, 2.68345455, 10.00000000, 0.00933207, 1.01007312),
    pars2 = c(100, 11, 1, 0, NA),
    datalist = Galapagos_datalist2,
    methode = "ode45",
    CS_version = list(model = 1, function_to_optimize = 'DAISIE'),
    abstolint = 1e-16,
    reltolint = 1e-10)
  testthat::expect_equal(loglik_CS01, loglik_CS11, tol = 5E-6)
})

test_that("DAISIE_loglik_CS_choice produces equivalent output for ODEINT RKCK54
          and deSolve lsodes", {
  #skip_if(Sys.getenv("CI") == "", message = "Run only on CI")
  #skip_on_cran()
  pars1 <- c(2.000, 2.700, 20.000, 0.009, 1.010)
  pars2 <- c(1.0e+02, 1.1e+01, 0.0e+00, 0.0e+00, NA)
  brts <- c(4.0000, 3.0282, 1.3227, 0.8223, 0.4286, 0.3462, 0.2450, 0.0808,
            0.0527, 0.0327, 0.0221, 0.1180, 0.0756, 0.0525, 0.0322, 0.0118)
  stac <- 2
  missnumspec <- 0
  loglik1 <- DAISIE_loglik_CS_choice(
    pars1 = pars1,
    pars2 = pars2,
    brts = brts,
    stac = stac,
    missnumspec = missnumspec
  )
  loglik2 <- DAISIE_loglik_CS_choice(
    pars1 = pars1,
    pars2 = pars2,
    brts = brts,
    stac = stac,
    missnumspec = missnumspec,
    methode = "odeint::runge_kutta_cash_karp54"
  )
  testthat::expect_equal(expected = loglik1, object = loglik2)
})

test_that("DAISIE_loglik_CS_choice produces equivalent
          output for ODEINT RKF78 and deSolve lsodes", {
  pars1 <- c(2.000, 2.700, 20.000, 0.009, 1.010)
  pars2 <- c(1.0e+02, 1.1e+01, 0.0e+00, 0.0e+00, NA)
  brts <- c(4.0000, 3.0282, 1.3227, 0.8223, 0.4286, 0.3462, 0.2450, 0.0808,
            0.0527, 0.0327, 0.0221, 0.1180, 0.0756, 0.0525, 0.0322, 0.0118)
  stac <- 2
  missnumspec <- 0
  CS_version <- 0
  # deSolve lsodes

  loglik1 <- expect_silent(
    DAISIE:::DAISIE_loglik_CS_choice(
      pars1 = pars1,
      pars2 = pars2,
      brts = brts,
      stac = stac,
      missnumspec = missnumspec,
      CS_version = CS_version
    )
  )
  # odeint RKF78
  loglik2 <- expect_silent(
    DAISIE:::DAISIE_loglik_CS_choice(
      pars1 = pars1,
      pars2 = pars2,
      brts = brts,
      stac = stac,
      missnumspec = missnumspec,
      CS_version = CS_version,
      methode = "odeint::runge_kutta_fehlberg78"
    )
  )
  testthat::expect_equal(expected = loglik1, object = loglik2)
})

