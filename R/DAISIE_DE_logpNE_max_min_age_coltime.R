#' @name DAISIE_DE_logpNE_max_min_age_coltime
#' @title Function to calculate the likelihood of observing a non-endemic lineage on the island
#' with minimum and maximum ages of colonization
#' @description This function calculates the log-likelihood of observing a non-endemic lineage on an island
#' for which the exact colonization time is unknown, but the maximum and minimum ages of colonization are known.
#' @inheritParams default_params_doc
#' @return The output is a numeric value representing the log-likelihood of observing a non-endemic singleton lineage
#' for which the minimum and maximum ages of colonization are given.
#' \item{logL1b}{ The log-likelihood value computed based on the differential equation system.}
#'
#' @export DAISIE_DE_logpNE_max_min_age_coltime

DAISIE_DE_logpNE_max_min_age_coltime <- function(brts,
                                                 pars1,
                                                 methode,
                                                 reltolint,
                                                 abstolint) {
  t0 <- brts[1]
  t1 <- brts[2]
  t2 <- brts[3]
  tp <- 0
  parameters <- pars1

  interval1 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dDm2 <- -(pars1[5] + pars1[1] + pars1[3] + pars1[4]) * Dm2
      dE <-  pars1[2] - (pars1[1] + pars1[2]) * E + pars1[1] * E^2
      list(c(dDm2, dE))
    })
  }

  # Time sequence for interval [t2, tp]
  time1 <- c(tp, t2)

  # Initial conditions
  initial_conditions1 <- c(Dm2 = 1, E = 0)

  # Solve the system for interval [t2, tp]
  solution1 <- deSolve::ode(y = initial_conditions1,
                            times = time1,
                            func = interval1,
                            parms = parameters,
                            method = methode,
                            rtol = reltolint,
                            atol = abstolint)

  interval2 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dDA <-  -pars1[4] * DA + pars1[4] * Dm2
      dDm1 <- -(pars1[5] + pars1[1] + pars1[3] + pars1[4]) * Dm1 +
        (pars1[5] * E + pars1[1] * E^2 + pars1[3]) * DA + pars1[4] * Dm2
      dDm2 <- -(pars1[5] + pars1[1] + pars1[3]) * Dm2 +
        (pars1[5] * E + pars1[1] * E^2 + pars1[3]) * DA
      dE <-  pars1[2] - (pars1[1] + pars1[2]) * E + pars1[1] * E^2
      list(c(dDA, dDm1, dDm2, dE))
    })
  }

  # Define system of equations for interval [t0, t1]
  interval3 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dDA <- -pars1[4] * DA + pars1[4] * Dm1
      dDm1 <- -(pars1[5] + pars1[1] + pars1[3]) * Dm1 + (pars1[5] * E + pars1[1] * E^2 + pars1[3]) * DA
      dE <- pars1[2] - (pars1[1] + pars1[2]) * E + pars1[1] * E^2
      list(c(dDA, dDm1, dE))
    })
  }

  # Initial conditions
  initial_conditions2 <- c(DA = 0, Dm1 = 0, Dm2 = solution1[, "Dm2"][[2]], E = solution1[, "E"][[2]])

  # Time sequence for interval [t1, t2]
  time2 <- c(t2, t1)

  # Solve the system for interval [t1, tp]
  solution2 <- deSolve::ode(y = initial_conditions2,
                            times = time2,
                            func = interval2,
                            parms = parameters,
                            method = methode,
                            rtol = reltolint,
                            atol = abstolint)

  # Initial conditions
  initial_conditions3 <- c(DA = solution2[, "DA"][[2]],
                           Dm1 = solution2[, "Dm1"][[2]],
                           E = solution2[, "E"][[2]])

  # Time sequence for interval [t0, t1]
  time3 <- c(t1, t0)

  # Solve the system for interval [t0, t1]
  solution3 <- deSolve::ode(y = initial_conditions3,
                            times = time3,
                            func = interval3,
                            parms = parameters,
                            method = methode,
                            rtol = reltolint,
                            atol = abstolint)

  # Extract log-likelihood
  L1 <- solution3[, "DA"][[2]]
  logL1b <- log(L1)
  return(logL1b)
}
