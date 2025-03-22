###############################################################################
### function to calculate the likelihood of observing an endemic non-singleton lineage at time t1
###############################################################################
### Using D-E approach


# pars1[1] corresponds to the Cladogenesis rate
# pars1[2] corresponds to the Extinction rate of endemic lineages
# pars1[3] corresponds to the Extinction rate of non-endemic lineages
# pars1[4] = corresponds to the Colonization rate
# pars1[5] = corresponds to the Anagenesis rate

DAISIE_DE_logpEC <- function(datalist,
                             i,
                             pars1,
                             methode,
                             reltolint,
                             abstolint) {

  t0 <- datalist[[1]]$island_age
  t1 <- datalist[[i]]$branching_times[2]
  t2 <- datalist[[i]]$branching_times[3]
  tp <- 0
  ti <- sort(datalist[[i]]$branching_times)
  ti <- ti[1:(length(ti)-2)]

  # Define system of equations for interval [t2, tp]
  interval1 <- function(t, state, pars1) {
    with(as.list(c(state, pars1)), {
      dD1 <- -(pars1[1] + pars1[2]) * D1 + 2 * pars1[1] * D1 * E1
      dD0 <- -pars1[4] * D0 + pars1[4] * Dm
      dDm <- -(pars1[5] + pars1[1] + pars1[3]) * Dm + (pars1[5] * E1 + pars1[1] * E1^2 + pars1[3]) * D0
      dE1 <- pars1[2] - (pars1[1] + pars1[2]) * E1 + pars1[1] * E1^2
      list(c(dD1, dD0, dDm, dE1))
    })
  }

  # Define system of equations for interval [t1, t2]
  interval2 <- function(t, state, pars1) {
    with(as.list(c(state, pars1)), {
      dD1 <- -(pars1[1] + pars1[2]) * D1 + 2 * pars1[1] * D1 * E1
      dD0 <- -pars1[4] * D0 + pars1[4] * Dm
      dDm <- -(pars1[5] + pars1[1] + pars1[3]) * Dm + (pars1[5] * E1 + pars1[1] * E1^2 + pars1[3]) * D0
      dDM <- -(pars1[5] + pars1[1] + pars1[3] + pars1[4]) * DM + (pars1[5] * D1 + 2 * pars1[1] * D1 * E1) * D0
      dE1 <- pars1[2] - (pars1[1] + pars1[2]) * E1 + pars1[1] * E1^2
      list(c(dD1, dD0, dDm, dDM, dE1))
    })
  }

  # Define system of equations for interval [t0, t1]
  interval3 <- function(t, state, pars1) {
    with(as.list(c(state, pars1)), {
      dD0 <- -pars1[4] * D0 + pars1[4] * Dm
      dDm <- -(pars1[5] + pars1[1] + pars1[3]) * Dm + (pars1[5] * E1 + pars1[1] * E1^2 + pars1[3]) * D0
      dE1 <- pars1[2] - (pars1[1] + pars1[2]) * E1 + pars1[1] * E1^2
      list(c(dD0, dDm, dE1))
    })
  }

  # Initial conditions
  number_of_species <- length(datalist[[i]]$branching_times) -1
  number_of_missing_species <- datalist[[i]]$missing_species
  ro <- number_of_species / (number_of_missing_species + number_of_species)
  initial_conditions1 <- c(D1 = ro, D0 = 1, Dm = 0, E1 = 1 - ro)

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
    initial_conditions1 <- c(D1 = pars1[1] * solution0[, "D1"][idx + 1] * solution1[, "D1"][2],
                             D0 = 1, Dm = 0, E1 = solution0[, "E1"][idx + 1])
  }

  # Initial conditions
  initial_conditions2 <- c(D1 = initial_conditions1["D1"][[1]],
                           D0 = solution0[, "D0"][length(ti) + 1],
                           Dm = solution0[, "Dm"][length(ti) + 1],
                           DM = initial_conditions1["D1"][[1]] * solution0[, "D0"][length(ti) + 1],
                           E1 = initial_conditions1["E1"][[1]])

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

  # Initial conditions
  initial_conditions3 <- c(D0 = pars1[4] * solution2[, "DM"][[2]],
                           Dm = pars1[4] * solution2[, "DM"][[2]],
                           E1 = solution2[, "E1"][[2]])

  # Time sequence for interval [t0, t1]
  time3 <- c(t1, t0)

  # Solve the system for interval [t0, t1]
  solution3 <- deSolve::ode(y = initial_conditions3,
                           times = time3,
                           func = interval3,
                           parms = pars1,
                           method = methode,
                           rtol = reltolint,
                           atol = abstolint)

  # Extract log-likelihood
  Lk <- solution3[, "D0"][[2]]
  logLkb <- log(Lk)
  return(logLkb)
}
