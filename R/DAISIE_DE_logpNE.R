###############################################################################
### function to calculate the likelihood of observing a non endemic lineages at time t1
###############################################################################
### Using D-E approach


# pars1[1] corresponds to the Cladogenesis rate
# pars1[2] corresponds to the Extinction rate of endemic lineages
# pars1[3] corresponds to the Extinction rate of non-endemic lineages
# pars1[4] = corresponds to the Colonization rate
# pars1[5] = corresponds to the Anagenesis rate

DAISIE_DE_logpNE <- function(datalist,i,
                             pars1,
                             methode,
                             reltolint,
                             abstolint) {
  t0 <- datalist[[1]]$island_age
  t1 <- datalist[[i]]$branching_times[2]
  tp <- 0
  parameters <- pars1

  # Define system of equations for interval [t1, tp]
  interval1 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dDM <- -(pars1[5] + pars1[1] + pars1[3] + pars1[4]) * DM
      dE1 <- pars1[2] - (pars1[1] + pars1[2]) * E1 + pars1[1] * E1^2
      list(c(dDM, dE1))
    })
  }

  # Define system of equations for interval [t0, t1]
  interval2 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dD0 <- -pars1[4] * D0 + pars1[4] * Dm
      dDm <- -(pars1[5] + pars1[1] + pars1[3]) * Dm + (pars1[5] * E1 + pars1[1] * E1^2 + pars1[3]) * D0
      dE1 <- pars1[2] - (pars1[1] + pars1[2]) * E1 + pars1[1] * E1^2
      list(c(dD0, dDm, dE1))
    })
  }

  # Set initial conditions
  initial_conditions1 <- c(DM = 1, E1 = 0)

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

  # Set initial conditions
  initial_conditions2 <- c(D0 = pars1[4] * solution1[, "DM"][[2]],
                           Dm = pars1[4] * solution1[, "DM"][[2]],
                           E1 = solution1[, "E1"][[2]])

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
  LM <- solution2[, "D0"][[2]]
  logLMb <- log(LM)
  return(logLMb)
}
