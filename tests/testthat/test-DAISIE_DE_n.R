test_that("DAISIE_DE_n gives the same result as DAISIE", {
  skip_on_cran()
  pars1 <- c(0.5,0.1,Inf,0.01,0.1)
  pars2 <- c(100,0,0,1)
  brts <- c(10,5, 3, 2)
  stac <- 2
  missnumspec <- 8
  methode <- 'odeint::runge_kutta_cash_karp54'
  abstolint <- 1E-12
  reltolint <- 1E-10
  verbose <- 1
  CS_version <- list(model = 1, function_to_optimize = 'DAISIE', sampling = 'n')
  loglik1 <- DAISIE_loglik(
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
  loglik_fun <- function(pars1, brts, missnumspec, methode) {
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
    return(loglikelihood)
  }
  #methode <- 'ode45'
  loglikelihood1 <- loglik_fun(pars1, brts, missnumspec, methode)
  #print(sprintf('%0.16f ',loglikelihood1))
  testthat::expect_equal(loglik1, loglikelihood1, tol = 1E-5) #-11.22540681977405085945, -11.22535738310295094777 / -11.2253390802734039

  pars1 <- c(0.25,0.1,Inf,0.01,0.1)
  loglik2 <- DAISIE_loglik(
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
  loglikelihood2 <- loglik_fun(pars1, brts, missnumspec, methode)
  #print(sprintf('%0.16f ',loglikelihood2))
  testthat::expect_equal(loglik2, loglikelihood2, tol = 1E-5) #-13.51826715621074193052, -13.51850213765633235141 / -13.5182968721357373

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
  loglik3 <- DAISIE_loglik(
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
  loglikelihood3 <- loglik_fun(pars1, brts, missnumspec, methode)
  #print(sprintf('%0.16f ',loglikelihood3))
  testthat::expect_equal(loglik3, loglikelihood3, tol = 1E-5) #-7.53368410913773978166, -7.53369312265383328509 / -7.5336902619524775
})
