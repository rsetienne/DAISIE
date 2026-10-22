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

  data("Galapagos_datalist", package = "DAISIE")
  pars_de <- c(2.546591, 2.678781, 2.678781, 0.009326754, 1.008583)

  for (i in c(4, 5)) {
    brts <- Galapagos_datalist[[i]]$branching_times
    stac <- Galapagos_datalist[[i]]$stac
    cat("clade", i, "stac", stac, "tips", length(brts) - 1, "\n")
    for (m in 0:10) {
      r <- tryCatch(
        DAISIE:::DAISIE_DE_loglik(pars1 = pars_de,
                                  brts = brts,
                                  missnumspec = m,
                                  stac = stac,
                                  methode = "odeint::runge_kutta_cash_karp54",
                                  sampling = "n"),
        error = function(e) paste("ERROR:", conditionMessage(e)))
    }
  }
  testthat::expect_equal(r, -16.51424, tol = 1e-4)
})

test_that("DAISIE_DE_n falls back from Cauchy integral to central and then
          forward differentiation", {
            brts <- c(10, 5, 3, 2)
            pars1 <- c(2.546591, 2.678781, 2.678781, 0.009326754, 1.008583)
            missnumspec <- 3
            # Wrapper around the real likelihood that fails on request:
            # - the Cauchy integral evaluates the likelihood at complex missnumspec,
            # - central differences evaluate it at x < 0, i.e. at missnumspec < 0,
            # - forward differences only evaluate it at x >= 0.
            make_function <- function(fail_complex = FALSE,
                                      fail_negative = FALSE,
                                      fail_nonzero = FALSE) {
              function(brts, missnumspec, stac, pars1, methode, reltolint, abstolint) {
                if (fail_complex && Im(missnumspec) != 0) stop("complex not allowed")
                if (fail_negative && Re(missnumspec) < 0) stop("negative not allowed")
                if (fail_nonzero && missnumspec != 0) stop("nonzero not allowed")
                DAISIE_DE_logpEC(brts = brts,
                                 missnumspec = missnumspec,
                                 stac = stac,
                                 pars1 = pars1,
                                 methode = methode,
                                 reltolint = reltolint,
                                 abstolint = abstolint)
              }
            }
            run <- function(DAISIE_DE_function) {
              DAISIE_DE_n(DAISIE_DE_function = DAISIE_DE_function,
                          brts = brts,
                          missnumspec = missnumspec,
                          stac = 2,
                          pars1 = pars1,
                          methode = "odeint::runge_kutta_cash_karp54",
                          reltolint = 1e-12,
                          abstolint = 1e-12)
            }

            # 1. Cauchy integral works: no messages
            expect_no_message(res_cauchy <- run(make_function()))

            # 2. Cauchy integral fails: central differentiation is used
            expect_message(res_central <- run(make_function(fail_complex = TRUE)),
                           "Cauchy integral failed")
            expect_equal(res_central, res_cauchy, tolerance = 1e-5)

            # 3. Cauchy integral and central differentiation fail: forward is used
            msgs <- testthat::capture_messages(
              res_forward <- run(make_function(fail_complex = TRUE,
                                               fail_negative = TRUE)))
            expect_length(msgs, 2)
            expect_match(msgs[1], "Cauchy integral failed")
            expect_match(msgs[2], "Central differentiation failed")
            expect_equal(res_forward, res_cauchy, tolerance = 1e-3)

            # 4. Everything fails: error
            suppressMessages(
              expect_error(run(make_function(fail_complex = TRUE,
                                             fail_negative = TRUE,
                                             fail_nonzero = TRUE)),
                           "all failed")
            )
          })

