test_that("rosenbrock4 and adams bashforth moulton1 work", {
  skip("Too slow to run")
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
  methode <- 'odeint::rosenbrock4'
  loglik_rb <- DAISIE_loglik_all(pars1, pars2, Galapagos_datalist_2types, methode = methode)
  methode <- 'odeint::adams_bashforth_moulton_1'
  DAISIE_CS_max_steps(100000000)
  DAISIE_abm_factor(0.000001)
  loglik_abm <- DAISIE_loglik_all(pars1, pars2, Galapagos_datalist_2types, methode = methode)
  testthat::expect_equal(loglik_lsodes, loglik_rb)
  testthat::expect_equal(loglik_lsodes, loglik_abm, tol = 1E-6)
})
