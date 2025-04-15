#' @name DAISIE_DE_logpEC
#' @title Function to calculate the likelihood of observing an endemic lineage
#' with the colonization time at t1.
#' @description This function calculates the log-likelihood of observing an endemic lineage
#' with the colonization time at t1.
#'
#' @inheritParams default_params_doc_DAISIE_DE
#' @return The output is a numeric value representing the log-likelihood of observing an endemic lineage
#' with the colonization time at t1.
#' \item{logLkb}{ The log-likelihood value computed based on a system of differential equations. }
#'
#' @examples
#'
#' # Select a dataset from a DAISIE package
#'
#' data(Galapagos_datalist)
#' datalist <- Galapagos_datalist
#'
#' # Select an endemic lineage in the dataset
#' i <- 4
#' # Define example parameters
#' pars1 <- c(0.2, 0.1, 0.05, 0.02, 0.03)
#'
#' # choose the method to solve the system of differential equations
#' log_likelihood <- DAISIE_DE_logpEC(datalist, i, pars1, methode = "lsodes", reltolint = 1e-16, abstolint = 1e-16)
#'
#' print(log_likelihood)
#'
#' @export DAISIE_DE_logpEC



DAISIE_DE_logpEC <- function(datalist,
                             i,
                             pars1,
                             methode,
                             reltolint,
                             abstolint) {

  brts = datalist[[i]]$branching_times
  missnumspec = datalist[[i]]$missing_species

  t0 <- brts[1]
  t1 <- brts[2]
  t2 <- brts[3]
  tp <- 0
  ti <- sort(brts)
  ti <- ti[1:(length(ti)-2)]

  # Define system of equations for interval [t2, tp]
  interval1 <- function(t, state, pars1) {
    with(as.list(c(state, pars1)), {
      dDE <- -(pars1[1] + pars1[2]) * DE + 2 * pars1[1] * DE * E
      dDA3 <- -pars1[4] * DA3 + pars1[4] * Dm3
      dDm3 <- -(pars1[5] + pars1[1] + pars1[3]) * Dm3 + (pars1[5] * E + pars1[1] * E^2 + pars1[3]) * DA3
      dE <- pars1[2] - (pars1[1] + pars1[2]) * E + pars1[1] * E^2
      list(c(dDE, dDA3, dDm3, dE))
    })
  }

  # Define system of equations for interval [t1, t2]
  interval2 <- function(t, state, pars1) {
    with(as.list(c(state, pars1)), {
      dDE <- -(pars1[1] + pars1[2]) * DE + 2 * pars1[1] * DE * E
      dDA3 <- -pars1[4] * DA3 + pars1[4] * Dm3
      dDm3 <- -(pars1[5] + pars1[1] + pars1[3]) * Dm3 + (pars1[5] * E + pars1[1] * E^2 + pars1[3]) * DA3
      dDm2 <- -(pars1[5] + pars1[1] + pars1[3] + pars1[4]) * Dm2 + (pars1[5] * DE + 2 * pars1[1] * DE * E) * DA3
      dE <- pars1[2] - (pars1[1] + pars1[2]) * E + pars1[1] * E^2
      list(c(dDE, dDA3, dDm3, dDm2, dE))
    })
  }

  # Define system of equations for interval [t0, t1]
  interval3 <- function(t, state, pars1) {
    with(as.list(c(state, pars1)), {
      dDA1 <- -pars1[4] * DA1 + pars1[4] * Dm1
      dDm1 <- -(pars1[5] + pars1[1] + pars1[3]) * Dm1 + (pars1[5] * E + pars1[1] * E^2 + pars1[3]) * DA1
      dE <- pars1[2] - (pars1[1] + pars1[2]) * E + pars1[1] * E^2
      list(c(dDA1, dDm1, dE))
    })
  }

  # Initial conditions
  number_of_species <- length(brts) -1
  number_of_missing_species <- missnumspec
  ro <- number_of_species / (number_of_missing_species + number_of_species)
  initial_conditions1 <- c(DE = ro, DA3 = 1, Dm3 = 0, E = 1 - ro)

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
                           DA3 = solution0[, "DA3"][length(ti) + 1],
                           Dm3 = solution0[, "Dm3"][length(ti) + 1],
                           Dm2 = initial_conditions1["DE"][[1]] * solution0[, "DA3"][length(ti) + 1],
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

  # Initial conditions
  initial_conditions3 <- c(DA1 = pars1[4] * solution2[, "Dm2"][[2]],
                           Dm1 = pars1[4] * solution2[, "Dm2"][[2]],
                           E = solution2[, "E"][[2]])

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
  Lk <- solution3[, "DA1"][[2]]
  logLkb <- log(Lk)
  return(logLkb)
}


