test_that("DAISIE_DE_n gives the same result as DAISIE", {
  skip_on_cran()
  pars1 <- c(0.5,0.1,Inf,0.01,0.1)
  pars2 <- c(200,0,0,1)
  brts1 <- c(10, 5, 3, 2)
  brts2 <- c(18, 7, 6)
  brts3 <- c(18, 7, 6, 5, 4, 3, 2, 1)
  stac <- 2
  missnumspec <- 18
  methode <- 'odeint::runge_kutta_cash_karp54'
  abstolint <- 1E-12
  reltolint <- 1E-10
  verbose <- 1
  CS_version <- list(model = 1, function_to_optimize = 'DAISIE', sampling = 'n')
  loglik1 <- DAISIE_loglik(
    pars1 = pars1,
    pars2 = pars2,
    brts = brts1,
    stac = stac,
    missnumspec = missnumspec,
    methode = methode,
    abstolint = abstolint,
    reltolint = reltolint,
    verbose = verbose,
    CS_version = CS_version
  )
  loglik2 <- DAISIE_loglik(
    pars1 = pars1,
    pars2 = pars2,
    brts = brts2,
    stac = stac,
    missnumspec = missnumspec,
    methode = methode,
    abstolint = abstolint,
    reltolint = reltolint,
    verbose = verbose,
    CS_version = CS_version
  )
  loglik3 <- DAISIE_loglik(
    pars1 = pars1,
    pars2 = pars2,
    brts = brts3,
    stac = stac,
    missnumspec = missnumspec,
    methode = methode,
    abstolint = abstolint,
    reltolint = reltolint,
    verbose = verbose,
    CS_version = CS_version
  )

  pars1[3] <- pars1[2]
  loglikelihood1 <- DAISIE_DE_loglik(pars1 = pars1,
                                     brts = brts1,
                                     missnumspec = missnumspec,
                                     stac = stac,
                                     methode = methode,                                     reltolint = 1e-15,
                                     abstolint = 1e-15,
                                     sampling = 'n')
  loglikelihood2 <- DAISIE_DE_loglik(pars1 = pars1,
                                     brts = brts2,
                                     missnumspec = missnumspec,
                                     stac = stac,
                                     methode = methode,
                                     reltolint = 1e-15,
                                     abstolint = 1e-15,
                                     sampling = 'n')
  loglikelihood3 <- DAISIE_DE_loglik(pars1 = pars1,
                                     brts = brts3,
                                     missnumspec = missnumspec,
                                     stac = stac,
                                     methode = methode,
                                     reltolint = 1e-15,
                                     abstolint = 1e-15,
                                     sampling = 'n')
  testthat::expect_equal(loglik1, loglikelihood1, tol = 1E-5)
  testthat::expect_equal(loglik2, loglikelihood2, tol = 1E-5)
  testthat::expect_equal(loglik3, loglikelihood3, tol = 1E-5)

  pars1 <- c(0.25,0.1,Inf,0.01,0.1)
  loglik4 <- DAISIE_loglik(
    pars1 = pars1,
    pars2 = pars2,
    brts = brts1,
    stac = stac,
    missnumspec = missnumspec,
    methode = methode,
    abstolint = abstolint,
    reltolint = reltolint,
    verbose = verbose,
    CS_version = CS_version
  )
  pars1[3] <- pars1[2]
  loglikelihood4 <- DAISIE_DE_loglik(pars1 = pars1,
                                     brts = brts1,
                                     missnumspec = missnumspec,
                                     stac = stac,
                                     methode = methode,
                                     reltolint = 1e-15,
                                     abstolint = 1e-15,
                                     sampling = 'n')
  testthat::expect_equal(loglik4, loglikelihood4, tol = 1E-5)

  pars1 <- c(0.5,0.1,Inf,0.01,0.1)
  pars2 <- c(100,0,0,1)
  brts <- c(10, 5)
  stac <- 2
  missnumspec <- 4
  methode <- 'odeint::runge_kutta_cash_karp54'
  abstolint <- 1E-12
  reltolint <- 1E-10
  verbose <- 1
  CS_version <- list(model = 1, function_to_optimize = 'DAISIE', sampling = 'n')
  loglik5 <- DAISIE_loglik(
    pars1 = pars1,
    pars2 = pars2,
    brts = brts,
    stac = stac,
    missnumspec = missnumspec,
    methode = methode,
    abstolint = abstolint,
    reltolint = reltolint,
    verbose = verbose,
    CS_version = CS_version
  )

  pars1[3] <- pars1[2]
  loglikelihood5 <- DAISIE_DE_loglik(pars1 = pars1,
                                     brts = brts,
                                     missnumspec = missnumspec,
                                     stac = stac,
                                     methode = methode,
                                     reltolint = 1e-15,
                                     abstolint = 1e-15,
                                     sampling = 'n')
  testthat::expect_equal(loglik5, loglikelihood5, tol = 1E-5)
})
