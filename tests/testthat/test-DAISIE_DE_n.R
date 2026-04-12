test_that("DAISIE_DE_n gives the same result as DAISIE", {
  pars1 <- c(0.5,0.1,Inf,0.01,0.1)
  pars2 <- c(100,0,0,1)
  brts <- c(10,5, 3, 2)
  stac <- 2
  missnumspec <- 2
  N_cheb <- 100
  methode <- 'odeint::runge_kutta_cash_karp54'
  abstolint <- 1E-12
  reltolint <- 1E-10
  verbose <- 1
  CS_version <- list(model = 1, function_to_optimize = 'DAISIE')
  loglik <- DAISIE_loglik(
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
  methode <- 'odeint::runge_kutta_cash_karp54'
  if (stac == 1 || stac == 4 || stac == 8) {
    loglikelihood <- DAISIE_DE_n(DAISIE_DE_function = DAISIE_DE_logpNE,
                                 brts = brts,
                                 pars1 = pars1,
                                 stac = stac,
                                 methode = methode,
                                 reltolint = 1e-15,
                                 abstolint = 1e-15)
  } else if (stac == 2 && length(brts) == 2 || stac == 3 && length(brts) == 2 || stac == 5 && length(brts) == 2 || stac == 9) {
    loglikelihood <- DAISIE_DE_n(DAISIE_DE_function = DAISIE_DE_logpES,
                                 brts = brts,
                                 missnumspec = missnumspec,
                                 stac = stac,
                                 pars1 = pars1,
                                 methode = methode,
                                 reltolint = 1e-15,
                                 abstolint = 1e-15)
  } else if (stac == 2 && length(brts) > 2 || stac == 3 && length(brts) > 2 || stac == 6) {
    loglikelihood <- DAISIE_DE_n(DAISIE_DE_function = DAISIE_DE_logpEC,
                                 brts = brts,
                                 missnumspec = missnumspec,
                                 stac = stac,
                                 pars1 = pars1,
                                 methode = methode,
                                 reltolint = 1e-15,
                                 abstolint = 1e-15)
  } else {
    stop("Unknown stac value: ", stac)
  }
  #print(loglik)
  #print(loglikelihood)
  testthat::expect_equal(loglik, loglikelihood, tol = 1E-5)
})
