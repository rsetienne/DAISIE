test_that("loglik IW various solver options give similar results", {
  # Test is not included in coverage due to issue with running loglik_IW
  # code from covr::package_coverage()
  testthat::skip_on_covr("Fails on covr")

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
