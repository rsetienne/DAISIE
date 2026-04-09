test_that("DAISIE_DE_n gives the same result as DAISIE", {
  pars1 <- c(0.5,0.1,Inf,0.01,0.1)
  pars2 <- c(100,0,0,1)
  brts <- c(10,5, 3, 2)
  stac <- 2
  missnumspec <- 3
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
  methode <- 'ode45'
  if (stac == 1) {
    loglikelihood <- DAISIE_DE_logpNE_max_age_coltime(brts,pars1,methode,reltolint,abstolint)
  } else if (stac == 2) {
    if (length(brts) == 2)
      loglikelihood <- DAISIE_DE_n(DAISIE_DE_logpES,brts,missnumspec,pars1,methode,reltolint,abstolint,N_cheb)
    else
      loglikelihood <- DAISIE_DE_n(DAISIE_DE_logpEC,brts,missnumspec,pars1,methode,reltolint,abstolint,N_cheb)

  } else if (stac == 3) {
    if (length(brts) == 2)
      loglikelihood <- DAISIE_DE_n(DAISIE_DE_logpES_mainland,brts,missnumspec,pars1,methode,reltolint,abstolint,N_cheb)
    else
      loglikelihood <- DAISIE_DE_n(DAISIE_DE_logpEC_mainland,brts,missnumspec,pars1,methode,reltolint,abstolint,N_cheb)
  } else if (stac == 4) {
    loglikelihood <- DAISIE_DE_logpNE(brts,pars1,methode,reltolint,abstolint)
  } else if (stac == 5) {
    loglikelihood <- DAISIE_DE_n(DAISIE_DE_logpES_max_age_coltime,brts,missnumspec,pars1,methode,reltolint,abstolint,N_cheb)
  } else if (stac == 6) {
    loglikelihood <- DAISIE_DE_n(DAISIE_DE_logpEC_max_age_coltime,brts,missnumspec,pars1,methode,reltolint,abstolint,N_cheb)
  } else if (stac == 7) {
    if (length(brts) == 2)
      loglikelihood <- DAISIE_DE_n(DAISIE_DE_logpES_max_age_coltime_and_mainland,brts,missnumspec,pars1,methode,reltolint,abstolint,N_cheb)
    else
      loglikelihood <- DAISIE_DE_n(DAISIE_DE_logpEC_max_age_coltime_and_mainland,brts,missnumspec,pars1,methode,reltolint,abstolint,N_cheb)
  } else if (stac == 8) {
    loglikelihood <- DAISIE_DE_logpNE_max_min_age_coltime(brts,pars1,methode,reltolint,abstolint)
  } else if (stac == 9) {
    loglikelihood <- DAISIE_DE_n(DAISIE_DE_logpES_max_min_age_coltime,brts,missnumspec,pars1,methode,reltolint,abstolint,N_cheb)
  } else {
    stop("Unknown stac value: ", stac)
  }
  #print(loglik)
  #print(loglikelihood)
  testthat::expect_equal(loglik, loglikelihood, tol = 1E-5)
})
