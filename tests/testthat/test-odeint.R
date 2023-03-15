test_that("odeint solvers give the same result as deSolve solvers", {
  utils::data(Galapagos_datalist_2types)
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
  pars2 <- c(40, 11, 0, 1)
  methode <- 'deSolve_R::lsodes'
  loglik_lsodes_R <- DAISIE_loglik_all(pars1, pars2, Galapagos_datalist_2types, methode = methode)
  methode <- 'deSolve_R::lsoda'
  loglik_lsoda_R <- DAISIE_loglik_all(pars1, pars2, Galapagos_datalist_2types, methode = methode)
  methode <- 'lsodes'
  loglik_lsodes <- DAISIE_loglik_all(pars1, pars2, Galapagos_datalist_2types, methode = methode)
  methode <- 'lsoda'
  loglik_lsoda <- DAISIE_loglik_all(pars1, pars2, Galapagos_datalist_2types, methode = methode)
  expect_equal(loglik_lsoda_R,loglik_lsoda)
  expect_equal(loglik_lsodes_R,loglik_lsodes)
  expect_equal(loglik_lsoda,loglik_lsodes)
  methode <- 'odeint::runge_kutta_cash_karp54'
  loglik_rkck54 <- DAISIE_loglik_all(pars1, pars2, Galapagos_datalist_2types, methode = methode)
  methode <- 'odeint::runge_kutta_fehlberg78'
  loglik_rkf78 <- DAISIE_loglik_all(pars1, pars2, Galapagos_datalist_2types, methode = methode)
  methode <- 'odeint::runge_kutta_dopri5'
  loglik_rkd5 <- DAISIE_loglik_all(pars1, pars2, Galapagos_datalist_2types, methode = methode)
  methode <- 'odeint::bulirsch_stoer'
  loglik_bs <- DAISIE_loglik_all(pars1, pars2, Galapagos_datalist_2types, methode = methode)
  expect_equal(loglik_lsodes, loglik_rkck54)
  expect_equal(loglik_lsodes, loglik_rkf78)
  expect_equal(loglik_lsodes, loglik_rkd5)
  expect_equal(loglik_lsodes, loglik_bs)


  pars1a <- pars1
  pars1a[6] <- Inf
  methode <- 'deSolve_R::lsoda'
  loglik_lsoda_R_Inf <- DAISIE_loglik_all(pars1a, pars2, Galapagos_datalist_2types, methode = methode)
  methode <- 'lsoda'
  loglik_lsoda_F_Inf <- DAISIE_loglik_all(pars1a, pars2, Galapagos_datalist_2types, methode = methode)
  expect_equal(loglik_lsoda_R_Inf,loglik_lsoda_F_Inf)

  methode <- 'deSolve_R::lsodes'
  loglik_lsodes_R_Inf <- DAISIE_loglik_all(pars1a, pars2, Galapagos_datalist_2types, methode = methode)
  methode <- 'lsodes'
  loglik_lsodes_F_Inf <- DAISIE_loglik_all(pars1a, pars2, Galapagos_datalist_2types, methode = methode)
  expect_equal(loglik_lsodes_R_Inf,loglik_lsodes_F_Inf)

  expect_equal(loglik_lsoda_R_Inf,loglik_lsoda_R, tol = 1E-4)
  expect_equal(loglik_lsodes_R_Inf,loglik_lsodes_R, tol = 1E-4)
})
