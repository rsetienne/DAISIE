test_that("IW and CS loglik is same when K = Inf", {
  skip_if(Sys.getenv("CI") == "", message = "Run only on CI")

  utils::data(Galapagos_datalist, package = "DAISIE")
  pars1 <- c(0.35, 0.3, Inf, 0.001, 0.3)
  pars2 <- c(120, 11, 0, 1)
  Galapagos_datalist_IW <- list()
  Galapagos_datalist_IW[[1]] <- Galapagos_datalist[[1]]
  Galapagos_datalist_IW[[1]]$not_present <- 1000
  Galapagos_datalist_IW[[2]] <- Galapagos_datalist[[2]]
  Galapagos_datalist_IW[[2]]$branching_times <- c(4, 2.9999999, 1.9998)
  Galapagos_datalist_IW[[2]]$stac <- 2
  Galapagos_datalist_IW[[3]] <- Galapagos_datalist[[3]]
  Galapagos_datalist_IW[[3]]$branching_times <- c(4, 1, 0.8)
  Galapagos_datalist_IW[[3]]$stac <- 2
  #Galapagos_datalist_IW <- Galapagos_datalist
  #for(i in 2:9) {
  #   Galapagos_datalist_IW[[i]]$branching_times <- c(4, 4 - 2*i*0.1,4 -2*i*0.1-0.1)
  #   Galapagos_datalist_IW[[i]]$stac <- 2
  #}

  #Galapagos_datalist_IW[[2]]$branching_times <- c(4, 3, 1.73)
  #Galapagos_datalist_IW[[2]]$stac <- 2
  #Galapagos_datalist_IW[[8]]$branching_times <- c(4, 2, 1.41)
  #Galapagos_datalist_IW[[8]]$stac <- 2

  Galapagos_datalist_IW <- DAISIE:::add_brt_table(Galapagos_datalist_IW)
  invisible(capture.output(
    loglik_IW <- DAISIE_loglik_IW(
      pars1 = pars1,
      pars2 = pars2,
      datalist = Galapagos_datalist_IW,
      methode = "odeint::runge_kutta_fehlberg78"
    )
  ))

  invisible(capture.output(
    loglik_CS <- DAISIE_loglik_CS(
      pars1 = pars1,
      pars2 = pars2,
      datalist = Galapagos_datalist_IW,
      methode = "odeint::runge_kutta_fehlberg78",
      CS_version = 1
    )
  ))
testthat::expect_equal(loglik_IW, loglik_CS, tol = 5E-6)
})

test_that("loglik IW various solver options give similar results", {
  skip_if(Sys.getenv("CI") == "", message = "Run only on CI")
  # Test is not included in coverage due to issue with running loglik_IW
  # code from covr::package_coverage()
  testthat::skip_on_covr()

  utils::data(frogs_datalist, package = "DAISIE")
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
  testthat::expect_equal(IW0,IW4, tolerance = 1E-4)
})
