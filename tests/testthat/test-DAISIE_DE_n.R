test_that("DAISIE_DE_n gives the same result as DAISIE", {
  skip_on_cran()
  pars1 <- c(0.5,0.1,Inf,0.01,0.1)
  pars2 <- c(200,0,0,1)
  brts <- list(c(10, 5, 3, 2),
               c(18, 7, 6),
               c(18, 7, 6, 5, 4, 3, 2, 1),
               c(18, 13, 12, 11, 10, 9, 8, 7, 6, 5))
  stac <- 2
  missnumspec <- 4
  methode <- 'odeint::runge_kutta_cash_karp54'
  abstolint <- 1E-12
  reltolint <- 1E-10
  verbose <- 1
  CS_version <- list(model = 1, function_to_optimize = 'DAISIE', sampling = 'n')
  loglik <- rep(0,length(brts))
  loglikelihood <- rep(0,length(brts))
  for(i in 1:length(brts)) {
    loglik[i] <- DAISIE_loglik(
      pars1 = pars1,
      pars2 = pars2,
      brts = brts[[i]],
      stac = stac,
      missnumspec = missnumspec,
      methode = methode,
      abstolint = abstolint,
      reltolint = reltolint,
      verbose = verbose,
      CS_version = CS_version
    )
    pars1new <- pars1
    pars1new[3] <- pars1new[2]
    loglikelihood[i] <- DAISIE_DE_loglik(pars1 = pars1new,
                                         brts = brts[[i]],
                                         missnumspec = missnumspec,
                                         stac = stac,
                                         methode = methode,                                     reltolint = 1e-15,
                                         abstolint = 1e-15,
                                         sampling = 'n')
    testthat::expect_equal(loglik[i], loglikelihood[i], tol = 1E-5)
  }

  pars1 <- c(0.25,0.1,Inf,0.01,0.1)
  for(i in 1:length(brts)) {
    loglik[i] <- DAISIE_loglik(
      pars1 = pars1,
      pars2 = pars2,
      brts = brts[[i]],
      stac = stac,
      missnumspec = missnumspec,
      methode = methode,
      abstolint = abstolint,
      reltolint = reltolint,
      verbose = verbose,
      CS_version = CS_version
    )
    pars1new <- pars1
    pars1new[3] <- pars1new[2]
    loglikelihood[i] <- DAISIE_DE_loglik(pars1 = pars1new,
                                         brts = brts[[i]],
                                         missnumspec = missnumspec,
                                         stac = stac,
                                         methode = methode,                                     reltolint = 1e-15,
                                         abstolint = 1e-15,
                                         sampling = 'n')
    testthat::expect_equal(loglik[i], loglikelihood[i], tol = 1E-5)
  }

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

  # data("Galapagos_datalist", package = "DAISIE")
  # pars_de <- c(2.546591, 2.678781, 2.678781, 0.009326754, 1.008583)
  # for (i in c(4, 5)) {
  #   brts <- Galapagos_datalist[[i]]$branching_times
  #   stac <- Galapagos_datalist[[i]]$stac
  #   cat("clade", i, "stac", stac, "tips", length(brts) - 1, "\n")
  #   for (m in 0:10) {
  #     r <- tryCatch(
  #       DAISIE:::DAISIE_DE_loglik(pars1 = pars_de,
  #                                 brts = brts,
  #                                 missnumspec = m,
  #                                 stac = stac,
  #                                 methode = "odeint::runge_kutta_cash_karp54",
  #                                 sampling = "n"),
  #       error = function(e) paste("ERROR:", conditionMessage(e)))
  #     cat("  m =", m, ":", format(r), "\n")
  #   }
  # }

})
