#' @name DAISIE_DE_logpEC_max_age_coltime_and_mainland
#' @title Function to calculate the likelihood of observing an endemic lineage on the island
#' with maximum age of colonization, and that coexists with its mainland ancestors. This is valid
#' for infinite K according to the DE equations.
#' @description This function calculates the log-likelihood of observing an endemic lineage on an island
#' with maximum age of colonization, and that coexists with its mainland ancestors. This is valid
#' for infinite K according to the DE equations.
#' @inheritParams default_params_doc
#' @return the loglikelihood
#' @noRd

DAISIE_DE_logpEC_max_age_coltime_and_mainland <- function(brts,
                                                          missnumspec,
                                                          pars1,
                                                          methode,
                                                          reltolint,
                                                          abstolint) {
  t0 <- brts[1]
  t1 <- brts[2]
  t2 <- brts[3]
  tp <- 0
  ti <- sort(brts)
  ti <- ti[1:(length(ti)-2)]
  parameters <- pars1


  # Define system of equations for interval [t2, tp]
  interval1 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dDE <- -(pars1[1] + pars1[2]) * DE + 2 * pars1[1] * DE * E
      dDA3 <- -pars1[4] * DA3 + pars1[4] * Dm3
      dDm3 <- -(pars1[5] + pars1[1] + pars1[3]) * Dm3 + (pars1[5] * E + pars1[1] * E^2 + pars1[3]) * DA3
      dE <- pars1[2] - (pars1[1] + pars1[2]) * E + pars1[1] * E^2
      list(c(dDE, dDA3, dDm3, dE))
    })
  }

  # Initial conditions
  number_of_species <- length(brts) - 1
  rho <- number_of_species / (missnumspec + number_of_species)

  initial_conditions1 <- c(DE = rho, DA3 = 0, Dm3 = 1, E = 1 - rho)

  # Define system of equations for interval [t1, t2]
  interval2 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dDE <- -(pars1[1] + pars1[2]) * DE + 2 * pars1[1] * DE * E
      dDA2 <- -pars1[4] * DA2 + pars1[4] * Dm2
      dDA3 <- -pars1[4] * DA3 + pars1[4] * Dm3
      dDm1 <- -(pars1[5] + pars1[1] + pars1[3] + pars1[4]) * Dm1 +
        (pars1[3] + pars1[5] * E + pars1[1] * E^2)* DA2 + pars1[4] * Dm2
      dDm2 <- -(pars1[5] + pars1[1] + pars1[3]) * Dm2 + (pars1[3] + pars1[5] * E + pars1[1] * E^2)* DA2 +
        (pars1[5] * DE + 2 * pars1[1] * DE * E ) * DA3
      dDm3 <- -(pars1[5] + pars1[1] + pars1[3]) * Dm3 + (pars1[3] + pars1[5] * E + pars1[1] * E^2) * DA3
      dE <- pars1[2] - (pars1[1] + pars1[2]) * E + pars1[1] * E^2
      list(c(dDE, dDA2, dDA3, dDm1, dDm2, dDm3, dE))
    })
  }

  interval3 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dDA1 <- -pars1[4] * DA1 + pars1[4] * Dm1
      dDm1 <- -(pars1[5] + pars1[1] + pars1[3]) * Dm1 + (pars1[5] * E + pars1[1] * E^2 + pars1[3]) * DA1
      dE <- pars1[2] - (pars1[1] + pars1[2]) * E + pars1[1] * E^2
      list(c(dDA1, dDm1, dE))
    })
  }

  solution0 <- deSolve::ode(y = initial_conditions1,
                            times = c(0, ti),
                            func = interval1,
                            parms = pars1,
                            method = methode,
                            rtol = reltolint,
                            atol = abstolint)

  # Time sequences for interval [t2, tp]
  times <- rbind(c(0, ti[1:(length(ti) - 1)]), ti)

  for (idx in 1:length(ti)) {
    # Time sequence idx in interval [t2, tp]
    time1 <- times[, idx]

    # Solve the system for interval [t2, tp]
    solution1 <- deSolve::ode(y = initial_conditions1,
                              times = time1,
                              func = interval1,
                              parms = pars1,
                              method = methode,
                              rtol = reltolint,
                              atol = abstolint)

    # Initial conditions
    initial_conditions1 <- c(DE = pars1[1] * solution0[, "DE"][idx + 1] * solution1[, "DE"][2],
                             DA3 = 1, Dm3 = 0, E = solution0[, "E"][idx + 1])
  }

  # Initial conditions
  initial_conditions2 <- c(DE = initial_conditions1["DE"][[1]],
                           DA2 = 0,
                           DA3 = solution0[, "DA3"][length(ti) + 1],
                           Dm1 = 0,
                           Dm2 = initial_conditions1["DE"][[1]] * solution0[, "DA3"][length(ti)+1],
                           Dm3 = solution0[, "Dm3"][length(ti) + 1],
                           E = initial_conditions1["E"][[1]])

  # Time sequence for interval [t1, t2]
  time2 <- c(t2, t1)

  # Solve the system for interval [t2, tp]
  solution2 <- deSolve::ode(y = initial_conditions2,
                            times = time2,
                            func = interval2,
                            parms = pars1,
                            method = methode,
                            rtol = reltolint,
                            atol = abstolint)

  # Time sequence for interval [t1, tp]
  time3 <- c(t1, t0)

  # Initial conditions
  initial_conditions3 <- c(DA1 = solution2[, "DA2"][[2]],
                           Dm1 = solution2[, "Dm1"][[2]],
                           E = solution2[, "E"][[2]])

  # Solve the system for interval [t0, t1]
  solution3 <- deSolve::ode(y = initial_conditions3,
                            times = time3,
                            func = interval3,
                            parms = pars1,
                            method = methode,
                            rtol = reltolint,
                            atol = abstolint)

  # Extract log-likelihood
  Lk <- (solution3[, "DA1"][[2]])
  logLkb <- log(Lk)
  return(logLkb)
}
