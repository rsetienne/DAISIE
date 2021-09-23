context("DAISIE_loglik_CS")

test_that("DAISIE_loglik_CS_choice produces correct output for CS_version 1", {
  pars1 <- c(2.000, 2.700, 20.000, 0.009, 1.010)
  pars2 <- c(1.0e+02, 1.1e+01, 0.0e+00, 0.0e+00, NA, 0.0e+00, 1.0e-04,
             1.0e-05, 1.0e-07, 3.0e+03, 9.5e-01, 9.8e-01)
  brts <- c(4.0000, 3.0282, 1.3227, 0.8223, 0.4286, 0.3462, 0.2450, 0.0808,
            0.0527, 0.0327, 0.0221, 0.1180, 0.0756, 0.0525, 0.0322, 0.0118)
  stac <- 2
  missnumspec <- 0
  loglik <- DAISIE_loglik_CS_choice(pars1 = pars1,
                                    pars2 = pars2,
                                    brts = brts,
                                    stac = stac,
                                    missnumspec = missnumspec)

  expect_true(is.numeric(loglik))
  expect_equal(loglik, -17.6535269346579)

})

test_that("DAISIE_loglik_CS_choice produces correct output for relaxed-rate
          model (CS_version = 2)", {
            pars1 <- c(2.000, 2.700, 20.000, 0.009, 1.010)
            pars2 <- c(1.0e+02, 1.1e+01, 0.0e+00, 0.0e+00, NA, 0.0e+00, 1.0e-04,
                       1.0e-05, 1.0e-07, 3.0e+03, 9.5e-01, 9.8e-01)
            brts <- c(4.0000, 3.0282, 1.3227, 0.8223, 0.4286, 0.3462, 0.2450, 0.0808,
                      0.0527, 0.0327, 0.0221, 0.1180, 0.0756, 0.0525, 0.0322, 0.0118)
            stac <- 2
            missnumspec <- 0
            CS_version <- list(model = 2,
                               relaxed_par = "cladogenesis",
                               sd = 1)

            invisible(capture.output(loglik <- DAISIE_loglik_CS_choice(pars1 = pars1,
                                                                       pars2 = pars2,
                                                                       brts = brts,
                                                                       stac = stac,
                                                                       missnumspec = missnumspec,
                                                                       CS_version = CS_version)))
            expect_true(is.numeric(loglik))
            expect_equal(loglik, -9.550184206825)

          })

test_that("DAISIE_loglik_CS_choice produces same output for CS_version = 0 (with M = 1) and CS_version = 1 ", {
  pars1 <- c(2.000, 2.700, 20.000, 0.009, 1.010)
  pars2 <- c(100, 11, 0, 0, NA, 0.0e+00, 1.0e-04,
             1.0e-05, 1.0e-07, 3.0e+03, 9.5e-01, 9.8e-01)
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

  expect_equal(loglik0,loglik1)
})

test_that("DAISIE_loglik_all produces correct output for relaxed-rate model", {
  utils::data(Galapagos_datalist)
  invisible(capture.output(suppressWarnings(
    loglik <- DAISIE::DAISIE_loglik_all(
      pars1 = c(2.55068735, 2.68345455, 10.00000000, 0.00933207, 1.01007312),
      pars2 = c(100, 0, 0, 0, NA),
      datalist = Galapagos_datalist,
      methode = "lsodes",
      CS_version = list(model = 2,
                        relaxed_par = "cladogenesis",
                        sd = 1),
      abstolint = 1e-16,
      reltolint = 1e-10
    )
  )))
  expect_true(is.numeric(loglik))
  expect_equal(loglik, -77.50300644907)
})

test_that("DAISIE_loglik produces correct output", {
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

test_that("DAISIE_loglik_all produces same output for CS_version 0 and 1 with and without conditioning", {
  utils::data(Galapagos_datalist)
  Galapagos_datalist2 <- Galapagos_datalist
  for(i in 2:9) {
    Galapagos_datalist2[[i]]$branching_times <- c(4, 4 - 2*i*0.1,4 -2*i*0.1-0.1)
    Galapagos_datalist2[[i]]$stac <- 2
  }
  Galapagos_datalist2 <- DAISIE:::add_brt_table(Galapagos_datalist2)
  loglik_CS00 <- DAISIE::DAISIE_loglik_all(
    pars1 = c(2.55068735, 2.68345455, 10.00000000, 0.00933207, 1.01007312),
    pars2 = c(100, 11, 0, 0, NA),
    datalist = Galapagos_datalist2,
    methode = "odeint::runge_kutta_fehlberg78",
    CS_version = 0,
    abstolint = 1e-16,
    reltolint = 1e-10)
  loglik_CS10 <- DAISIE::DAISIE_loglik_all(
    pars1 = c(2.55068735, 2.68345455, 10.00000000, 0.00933207, 1.01007312),
    pars2 = c(100, 11, 0, 0, NA),
    datalist = Galapagos_datalist2,
    methode = "ode45",
    CS_version = 1,
    abstolint = 1e-16,
    reltolint = 1e-10)
  testthat::expect_equal(loglik_CS00, loglik_CS10, tol = 5E-6)
  loglik_CS01 <- DAISIE::DAISIE_loglik_all(
    pars1 = c(2.55068735, 2.68345455, 10.00000000, 0.00933207, 1.01007312),
    pars2 = c(100, 11, 1, 0, NA),
    datalist = Galapagos_datalist2,
    methode = "odeint::runge_kutta_fehlberg78",
    CS_version = 0,
    abstolint = 1e-16,
    reltolint = 1e-10)
  loglik_CS11 <- DAISIE::DAISIE_loglik_all(
    pars1 = c(2.55068735, 2.68345455, 10.00000000, 0.00933207, 1.01007312),
    pars2 = c(100, 11, 1, 0, NA),
    datalist = Galapagos_datalist2,
    methode = "ode45",
    CS_version = 1,
    abstolint = 1e-16,
    reltolint = 1e-10)
  testthat::expect_equal(loglik_CS01, loglik_CS11, tol = 5E-6)
})

test_that("DAISIE_loglik_CS_choice produces equivalent output for ODEINT RKCK54
          and deSolve lsodes", {
  pars1 <- c(2.000, 2.700, 20.000, 0.009, 1.010)
  pars2 <- c(1.0e+02, 1.1e+01, 0.0e+00, 0.0e+00, NA, 0.0e+00, 1.0e-04,
             1.0e-05, 1.0e-07, 3.0e+03, 9.5e-01, 9.8e-01)
  brts <- c(4.0000, 3.0282, 1.3227, 0.8223, 0.4286, 0.3462, 0.2450, 0.0808,
            0.0527, 0.0327, 0.0221, 0.1180, 0.0756, 0.0525, 0.0322, 0.0118)
  stac <- 2
  missnumspec <- 0
  loglik1 <- DAISIE:::DAISIE_loglik_CS_choice(pars1 = pars1,
                                              pars2 = pars2,
                                              brts = brts,
                                              stac = stac,
                                              missnumspec = missnumspec)
  loglik2 <- DAISIE:::DAISIE_loglik_CS_choice(pars1 = pars1,
                                              pars2 = pars2,
                                              brts = brts,
                                              stac = stac,
                                              missnumspec = missnumspec,
                                              methode = "odeint::runge_kutta_cash_karp54")
  expect_equal(expected = loglik1, object = loglik2)
})



test_that("DAISIE_loglik_CS_choice produces equivalent
          output for ODEINT RKF78 and deSolve lsodes", {
  pars1 <- c(2.000, 2.700, 20.000, 0.009, 1.010)
  pars2 <- c(1.0e+02, 1.1e+01, 0.0e+00, 0.0e+00, NA, 0.0e+00, 1.0e-04,
             1.0e-05, 1.0e-07, 3.0e+03, 9.5e-01, 9.8e-01)
  brts <- c(4.0000, 3.0282, 1.3227, 0.8223, 0.4286, 0.3462, 0.2450, 0.0808,
            0.0527, 0.0327, 0.0221, 0.1180, 0.0756, 0.0525, 0.0322, 0.0118)
  stac <- 2
  missnumspec <- 0
  CS_version <- 0
  # deSolve lsodes
  loglik1 <- DAISIE_loglik_CS_choice(pars1 = pars1,
                                     pars2 = pars2,
                                     brts = brts,
                                     stac = stac,
                                     missnumspec = missnumspec,
                                     CS_version = CS_version)
  # odeint RKF78
  loglik2 <- DAISIE_loglik_CS_choice(pars1 = pars1,
                                     pars2 = pars2,
                                     brts = brts,
                                     stac = stac,
                                     missnumspec = missnumspec,
                                     CS_version = CS_version,
                                     methode = "odeint::runge_kutta_fehlberg78")
  expect_equal(expected = loglik1, object = loglik2)
})

test_that("DAISIE_loglik_CS_choice produces equivalent output for ontogeny deSolve lsodes and odeint RKF78", {

  lac0 <- 2.000
  mu0 <- 2.700
  K0 <- 20.000
  gam0 <- 0.009
  laa0 <- 1.010
  d <- 0
  x <- 0
  area_pars <- c(
    max_area = 1,
    current_area = 1,
    proportional_peak_t = 0,
    total_island_age = 4,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0
  )
  island_ontogeny <- 0
  sea_level <- 0
  totaltime <- 4
  peak <- 1

  pars1_time_dep <- c(
    lac0,
    mu0,
    K0,
    gam0,
    laa0,
    d,
    x,
    area_pars,
    island_ontogeny,
    sea_level,
    totaltime,
    peak
  )
  pars2 <- c(1.0e+02, 1.1e+01, 0.0e+00, 0.0e+00, NA, 0.0e+00, 1.0e-04,
             1.0e-05, 1.0e-07, 3.0e+03, 9.5e-01, 9.8e-01)
  brts <- c(4.0000, 3.0282, 1.3227, 0.8223, 0.4286, 0.3462, 0.2450,
            0.0808, 0.0527, 0.0327, 0.0221, 0.1180, 0.0756, 0.0525,
            0.0322, 0.0118)


  pars1_const_rate <- c(2.000, 2.700, 20.000, 0.009, 1.010)

  stac <- 2
  missnumspec <- 0
  CS_version <- 0
  # deSolve lsodes time dep function with A = 1
  loglik1 <- DAISIE_loglik_CS_choice(pars1 = pars1_time_dep,
                                     pars2 = pars2,
                                     brts = brts,
                                     stac = stac,
                                     missnumspec = missnumspec,
                                     CS_version = CS_version)

  # odeint RKF78 constant rate function

  loglik2 <- DAISIE_loglik_CS_choice(pars1 = pars1_const_rate,
                                     pars2 = pars2,
                                     brts = brts,
                                     stac = stac,
                                     missnumspec = missnumspec,
                                     CS_version = CS_version,
                                     methode = "odeint::runge_kutta_fehlberg78")

  # deSolve lsodes constant rate function
  loglik3 <- DAISIE_loglik_CS_choice(pars1 = pars1_const_rate,
                                     pars2 = pars2,
                                     brts = brts,
                                     stac = stac,
                                     missnumspec = missnumspec,
                                     CS_version = CS_version)

  expect_equal(expected = loglik1, object = loglik2)
  expect_equal(expected = loglik1, object = loglik3)
})
