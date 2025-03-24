###############################################################################
### function to calculate the likelihood of observing a non-endemic lineage
### with the max and min age colonization time t1
###############################################################################
### Using D-E approach


# pars1[1] corresponds to the Cladogenesis rate
# pars1[2] corresponds to the Extinction rate of endemic lineages
# pars1[3] corresponds to the Extinction rate of non-endemic lineages
# pars1[4] = corresponds to the Colonization rate
# pars1[5] = corresponds to the Anagenesis rate

### Using D-E approach
DAISIE_DE_logpNE_max_min_age_coltime <- function(datalist,
                                             i,
                                             pars1,
                                             methode,
                                             reltolint,
                                             abstolint) {
  t0 <- datalist[[i]]$branching_times[1]
  t1 <- datalist[[i]]$branching_times[2]
  t2 <- datalist[[i]]$branching_times[3]
  tp <- 0
  parameters <- pars1

  interval1 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {

      dDM <- -(pars1[5] + pars1[1] + pars1[2] + pars1[4]) * DM

      dE1 <-  pars1[2] - (pars1[1] + pars1[2]) * E1 + pars1[1] * E1^2

      list(c(dDM, dE1))
    })
  }

  # Time sequence for interval [t2, tp]
  time1 <- c(tp, t2)

  # Initial conditions
  initial_conditions1 <- c(DM = 1, E1 = 0)

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

      dD0 <-  -pars1[4] * D0 + pars1[4] * DM

      dDm <- -(pars1[5] + pars1[1] + pars1[2] + pars1[4]) * Dm +
        (pars1[5] * E1 + pars1[1] * E1^2 + pars1[2]) * D0 + pars1[4] * DM

      dDM <- -(pars1[5] + pars1[1] + pars1[2]) * DM +
        (pars1[5] * E1 + pars1[1] * E1^2 + pars1[2]) * D0

      dE1 <-  pars1[2] - (pars1[1] + pars1[2]) * E1 + pars1[1] * E1^2

      list(c(dD0, dDm, dDM, dE1))
    })
  }



  # Define system of equations for interval [t0, t1]
  interval3 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dD0 <- -pars1[4] * D0 + pars1[4] * Dm

      dDm <- -(pars1[5] + pars1[1] + pars1[2]) * Dm + (pars1[5] * E1 + pars1[1] * E1^2 + pars1[2]) * D0

      dE1 <- pars1[2] - (pars1[1] + pars1[2]) * E1 + pars1[1] * E1^2

      list(c(dD0, dDm, dE1))
    })
  }

  # Initial conditions
  initial_conditions2 <- c(D0 = 0, Dm = 0, DM = solution1[, "DM"][[2]], E1 = solution1[, "E1"][[2]])

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
  initial_conditions3 <- c(D0 = solution2[, "D0"][[2]],
                           Dm = solution2[, "Dm"][[2]],
                           E1 = solution2[, "E1"][[2]])

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
  L1 <- solution3[, "D0"][[2]]
  logL1b <- log(L1)
  return(logL1b)
}


i <- 4
pars2 <- c(100, 11,0,2)
datalist[[i]]$branching_times <- c(4, 3, 1)

DAISIE_DE_logpNE_max_min_age_coltime(datalist,
                                     i,
                                     pars1,
                                     methode = "lsodes",
                                     reltolint = 1e-16,
                                     abstolint = 1e-16)

DAISIE:::DAISIE_loglik_CS_choice(pars1 = pars1,
                                 pars2 = pars2,
                                 datalist = datalist,
                                 brts = datalist[[i]]$branching_times,
                                 stac = 8,
                                 missnumspec = datalist[[i]]$missing_species,
                                 methode = "odeint::runge_kutta_fehlberg78")

