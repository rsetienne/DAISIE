#' @name DAISIE_DE_logpES_mainland
#' @title Function to calculate the likelihood of observing an endemic singleton
#' lineage on the island that coexists with its mainland ancestor. This is valid
#' for infinite K according to the DE equations.
#' @description This function calculates the log-likelihood of observing an
#' endemic singleton lineage on an island that coexists with its mainland
#' ancestor. This is valid for infinite K according to the DE equations.
#' @inheritParams default_params_doc
#' @return the loglikelihood
#' @noRd

DAISIE_DE_logpES_mainland <- function(brts,
                                      missnumspec,
                                      pars1,
                                      methode,
                                      reltolint,
                                      abstolint) {
  t0 <- brts[1]
  t1 <- brts[2]
  tp <- 0
  parameters <- pars1

  # Define system of equations for interval [t1, tp]
  interval1 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dDE <- -(pars1[1] + pars1[2]) * DE + 2 * pars1[1] * DE * E
      dDA3 <- -pars1[4] * DA3 + pars1[4] * Dm
      dDm3 <- -(pars1[5] + pars1[1] + pars1[3]) * Dm + (pars1[5] * E + pars1[1] * E^2 + pars1[3]) * DA3
      dDm2 <- -(pars1[5] + pars1[1] + pars1[3] + pars1[4]) * Dm2 + (pars1[5] * D1 + 2 * pars1[1] * DE * E) * DA3
      dE <- pars1[2] - (pars1[1] + pars1[2]) * E + pars1[1] * E^2
      list(c(dDE, dDA3, dDm3, dDm2, dE))
    })
  }

  # Define system of equations for interval [t0, t1]
  interval2 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dDA1 <- -pars1[4] * DA1 + pars1[4] * Dm1
      dDm1 <- -(pars1[5] + pars1[1] + pars1[3]) * Dm1 + (pars1[5] * E + pars1[1] * E^2 + pars1[3]) * DA1
      dE <- pars1[2] - (pars1[1] + pars1[2]) * E + pars1[1] * E^2
      list(c(dDA1, dDm1, dE))
    })
  }

  number_of_species <- length(brts) - 1
  rho <- number_of_species / (missnumspec + number_of_species)

  # Initial conditions
  initial_conditions1 <- c(DE = rho, DA3 = 0, Dm3 = 1, Dm2 = 0, E = 1 - rho)

    # Time sequence for interval [t1, tp]
  time1 <- c(tp, t1)

  # Solve the system for interval [t1, tp]
  solution1 <- deSolve::ode(y = initial_conditions1,
                            times = time1,
                            func = interval1,
                            parms = parameters,
                            method = methode,
                            rtol = reltolint,
                            atol = abstolint)

  # Initial conditions
  initial_conditions2 <- c(DA1 = pars1[4] * solution1[, "Dm2"][[2]],
                           Dm1 = pars1[4] * solution1[, "Dm2"][[2]],
                           E = solution1[, "E"][[2]])

  # Time sequence for interval [t0, t1]
  time2 <- c(t1, t0)

  # Solve the system for interval [t0, t1]
  solution2 <- deSolve::ode(y = initial_conditions2,
                            times = time2,
                            func = interval2,
                            parms = parameters,
                            method = methode,
                            rtol = reltolint,
                            atol = abstolint)

  # Extract log-likelihood
  L1 <- solution2[, "DA1"][[2]]
  logL1b <- log(L1)
  return(logL1b)
}
