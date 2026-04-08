test_that("DAISIE_logp0 is correct", {
  data("Galapagos_datalist", package = "DAISIE")
  datalist <- Galapagos_datalist

  parameter <- c(2.546591, 2.678781, Inf, 0.009326754, 1.008583)

  res1 <- DAISIE:::DAISIE_DE_logp0(
    island_age = datalist[[1]]$island_age,
    pars1 = c(2.546591, 2.678781, 2.678781, 0.009326754, 1.008583),
    methode = "lsodes",
    reltolint = 1e-15,
    abstolint = 1e-15)

  res2 <- DAISIE:::DAISIE_DE_logp0(island_age = datalist[[1]]$island_age,
                                   parameter,
                                   abstolint = 1e-15,
                                   reltolint = 1e-15,
                                   methode = "ode45"
  )

  testthat::expect_equal(res1, res2, tolerance = 1e-6)

  res3 <- DAISIE:::DAISIE_DE_logp0(island_age = datalist[[1]]$island_age,
                                   parameter,
                                   abstolint = 1e-15,
                                   reltolint = 1e-15,
                                   methode = "odeint::runge_kutta_cash_karp54"
  )
  testthat::expect_equal(res3, res2, tolerance = 1e-6)
})

test_that("logpEC", {
  data("Galapagos_datalist", package = "DAISIE")
  datalist <- Galapagos_datalist

  i <- 4
  brts <- datalist[[i]]$branching_times


  parameter <- c(2.546591, 2.678781, Inf, 0.009326754, 1.008583)

  res1 <- DAISIE:::DAISIE_DE_logpEC(brts,
                                    missnumspec = 0,
                                    stac        = 2,
                                    pars1       = parameter,
                                    abstolint   = 1e-15,
                                    reltolint   = 1e-15,
                                    methode     = "ode45")

  res2 <- DAISIE:::DAISIE_loglik_CS_choice(pars1 = c(2.546591, 2.678781, Inf, 0.009326754, 1.008583),
                                           pars2 = c(100, 11, 0, 2),
                                           brts = brts,
                                           stac = 2,
                                           missnumspec = 0,
                                           datalist = datalist)
  testthat::expect_equal(res1, res2)

  res3 <- DAISIE:::DAISIE_DE_logpEC(brts,
                                    missnumspec = 0,
                                    stac = 2,
                                    pars1 = parameter,
                                    abstolint  = 1e-15,
                                    reltolint  = 1e-15,
                                    methode = "odeint::runge_kutta_cash_karp54")
  testthat::expect_equal(res3, res1)
})

test_that("logpES", {

  data("Galapagos_datalist", package = "DAISIE")
  datalist <- Galapagos_datalist

  i <- 6
  brts <- datalist[[i]]$branching_times

  parameter <- c(2.546591, 2.678781, Inf, 0.009326754, 1.008583)

  res1 <- DAISIE:::DAISIE_DE_logpES(brts,
                                    missnumspec = 0,
                                    stac = 2,
                                    methode = "ode45",
                                    pars1 = parameter,
                                    abstolint  = 1e-15,
                                    reltolint  = 1e-15)

  res2 <- DAISIE:::DAISIE_loglik_CS_choice(pars1 = c(2.546591, 2.678781, Inf, 0.009326754, 1.008583),
                                           pars2 = c(100, 11, 0, 2),
                                           brts = brts,
                                           stac = 2,
                                           missnumspec = 0,
                                           datalist = datalist)

  testthat::expect_equal(res1, res2)

  res3 <- DAISIE:::DAISIE_DE_logpES(brts,
                                    missnumspec = 0,
                                    stac = 2,
                                    pars1 = parameter,
                                    abstolint  = 1e-15,
                                    reltolint  = 1e-15,
                                    methode = "odeint::runge_kutta_cash_karp54")
  testthat::expect_equal(res3, res1)

  ## stac 9
  data("Galapagos_datalist", package = "DAISIE")
  datalist <- Galapagos_datalist


  brts <- c(8, 5, 3)


  parameter <- c(2.546591, 2.678781, Inf, 0.009326754, 1.008583)

  res1 <- DAISIE:::DAISIE_DE_logpES(brts,
                                    missnumspec = 0,
                                    stac = 9,
                                    pars1 = parameter,
                                    abstolint  = 1e-15,
                                    reltolint  = 1e-15,
                                    methode = "ode45")

  res2 <- DAISIE:::DAISIE_loglik_CS_choice(pars1 = parameter,
                                           pars2 = c(100, 11, 0, 2),
                                           brts = brts,
                                           stac = 9,
                                           missnumspec = 0,
                                           datalist = datalist)
  testthat::expect_equal(res1, res2)

  res1 <- DAISIE:::DAISIE_DE_logpES(brts,
                                    missnumspec = 0,
                                    stac = 9,
                                    pars1 = parameter,
                                    abstolint  = 1e-15,
                                    reltolint  = 1e-15,
                                    methode = "odeint::runge_kutta_cash_karp54")
  testthat::expect_equal(res1, res2)
})

test_that("logpNE", {

  data("Galapagos_datalist", package = "DAISIE")
  datalist <- Galapagos_datalist

  i <- 3
  brts <- datalist[[i]]$branching_times


  parameter <- c(2.546591, 2.678781, Inf, 0.009326754, 1.008583)

  res1 <- DAISIE:::DAISIE_DE_logpNE(brts,
                                    pars1 = parameter,
                                    stac = 4,
                                    abstolint  = 1e-15,
                                    reltolint  = 1e-15,
                                    methode = "ode45")

  res2 <- DAISIE:::DAISIE_loglik_CS_choice(pars1 = parameter,
                                           pars2 = c(100, 11, 0, 2),
                                           brts = brts,
                                           stac = 4,
                                           missnumspec = 0,
                                           datalist = datalist)
  testthat::expect_equal(res1, res2)

  res2 <- DAISIE:::DAISIE_DE_logpNE(brts,
                                    pars1 = parameter,
                                    stac = 4,
                                    abstolint  = 1e-15,
                                    reltolint  = 1e-15,
                                    methode = "odeint::runge_kutta_cash_karp54")
  testthat::expect_equal(res1, res2)

  # stac 8
  data("Galapagos_datalist", package = "DAISIE")
  datalist <- Galapagos_datalist


  brts <- c(5, 4, 3)


  parameter <- c(2.546591, 2.678781, Inf, 0.009326754, 1.008583)

  res1 <- DAISIE:::DAISIE_DE_logpNE(brts,
                                    pars1 = parameter,
                                    stac = 8,
                                    abstolint  = 1e-15,
                                    reltolint  = 1e-15,
                                    methode = "ode45")

  res2 <- DAISIE:::DAISIE_loglik_CS_choice(pars1 = c(2.546591, 2.678781, Inf, 0.009326754, 1.008583),
                                           pars2 = c(100, 11, 0, 2),
                                           brts = brts,
                                           stac = 8,
                                           missnumspec = 0,
                                           datalist = datalist)
  testthat::expect_equal(res1, res2)

  res1 <- DAISIE:::DAISIE_DE_logpNE(brts,
                                    pars1 = parameter,
                                    stac = 8,
                                    abstolint  = 1e-15,
                                    reltolint  = 1e-15,
                                    methode = "odeint::runge_kutta_cash_karp54")
  testthat::expect_equal(res1, res2)
})
