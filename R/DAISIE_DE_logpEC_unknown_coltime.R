#' @name DAISIE_DE_logpEC_unknown_coltime
#' @title Function to calculate the likelihood of observing an endemic lineage
#' on the island with unknown colonization time. This is valid for infinite K
#' according to the DE equations.
#' @description This function calculates the log-likelihood of observing an
#' endemic lineage on an island for which the exact colonization time is unknown.
#' This is valid for infinite K according to the DE equations.
#' @inheritParams default_params_doc
#' @return the loglikelihood
#' @noRd

DAISIE_DE_logpEC_unknown_coltime <- function(brts,
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
      dD1 <- -(pars1[1] + pars1[2]) * D1 + 2 * pars1[1] * D1 * E1
      dD0 <- -pars1[4] * D0 + pars1[4] * Dm
      dDm <- -(pars1[5] + pars1[1] + pars1[2]) * Dm + (pars1[5] * E1 + pars1[1] * E1^2 + pars1[2]) * D0
      dE1 <- pars1[2] - (pars1[1] + pars1[2]) * E1 + pars1[1] * E1^2
      list(c(dD1, dD0, dDm, dE1))
    })
  }

  # Define system of equations for interval [t1, t2]
  interval2 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dD1 <- -(pars1[1] + pars1[2] ) * D1 + 2 * pars1[1] * D1 * E1
      dD0m <- -pars1[4] * D0m + pars1[4] * Dm
      dD0M <- -pars1[4] * D0M + pars1[4]*DM
      dDm <- -(pars1[5] + pars1[1] + pars1[2]) * Dm + (pars1[5] * E1 + pars1[1] * E1^2 + pars1[2]) * D0m
      dDM <- -(pars1[5] + pars1[1] + pars1[2]) * DM + (pars1[5] * E1 + pars1[1] * E1^2 + pars1[2])* D0M +
        (pars1[5] * D1 + 2 * pars1[1] * D1 * E1 ) * D0m
      dE1 <- pars1[2] - (pars1[1] + pars1[2]) * E1 + pars1[1] * E1^2
      list(c(dD1, dD0m, dD0M, dDm, dDM, dE1))
    })
  }

  # Initial conditions
  number_of_species <- length(brts) - 1
  rho <- number_of_species / (missnumspec + number_of_species)

  initial_conditions1 <- c(D1 = rho, D0 = 1, Dm = 0, E1 = 1 - rho)

  solution0 <- deSolve::ode(y = initial_conditions1,
                            times = c(0, ti),
                            func = interval1,
                            parms = parameters,
                            method = "lsodes",
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
                              parms = parameters,
                              method = "lsodes",
                              rtol = reltolint,
                              atol = abstolint)

    # Initial conditions
    initial_conditions1 <- c(D1 = pars1[1] * solution0[, "D1"][idx + 1] * solution1[, "D1"][2],
                             D0 = 1, Dm = 0, E1 = solution0[, "E1"][idx + 1])
  }

  # Initial conditions
  initial_conditions2 <- c(D1 = initial_conditions1["D1"][[1]],
                           D0m = solution0[, "D0"][length(ti) + 1],
                           D0M = 0,
                           Dm = solution0[, "Dm"][length(ti) + 1],
                           DM = initial_conditions1["D1"][[1]] * solution0[, "D0"][length(ti)+1],
                           E1 = initial_conditions1["E1"][[1]])

  # Time sequence for interval [t1, t2]
  time2 <- c(t2, t1)

  # Solve the system for interval [t2, tp]
  solution2 <- deSolve::ode(y = initial_conditions2,
                            times = time2,
                            func = interval2,
                            parms = parameters,
                            method = "lsodes",
                            rtol = reltolint,
                            atol = abstolint)


  # Extract log-likelihood
  Lk <- (solution2[, "D0M"][[2]])
  logLkb <- log(Lk)
  return(logLkb)
}
