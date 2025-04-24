#' @name DAISIE_DE_logpNE
#' @title Function to calculate the likelihood of observing a non-endemic lineage
#' with the colonization time at t1.
#' @description This function calculates the log-likelihood of observing a non-endemic lineage
#' with the colonization time at t1.
#'
#' @inheritParams default_params_doc_DAISIE_DE
#' @return The output is a numeric value representing the log-likelihood of observing a non-endemic lineage
#' with the colonization time at t1.
#' \item{logLMb}{ The log-likelihood value computed based on the differential equation system. }
#'
#' @examples
#'
#' # Select a dataset from a DAISIE package
#'
#' data(Galapagos_datalist)
#' datalist <- Galapagos_datalist
#'
#' # Select a non-endemic lineage in the dataset
#' i <- 3
#' # Define example parameters
#' pars1 <- c(0.2, 0.1, 0.05, 0.02, 0.03)
#'
#' # choose the method to solve the system of differential equations
#' log_likelihood <- DAISIE_DE_logpNE(datalist, i, pars1, methode = "lsodes", reltolint = 1e-16, abstolint = 1e-16)
#'
#' print(log_likelihood)
#'
#' @export DAISIE_DE_logpNE




DAISIE_DE_logpNE <- function(datalist,
                             i,
                             pars1,
                             methode,
                             reltolint,
                             abstolint) {

  brts = datalist[[i]]$branching_times
  missnumspec = datalist[[i]]$missing_species

  t0 <- brts[1]
  t1 <- brts[2]
  tp <- 0
  parameters <- pars1

  # Define system of equations for interval [t1, tp]
  interval1 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dDm2 <- -(pars1[5] + pars1[1] + pars1[3] + pars1[4]) * Dm2
      dE <- pars1[2] - (pars1[1] + pars1[2]) * E + pars1[1] * E^2
      list(c(dDm2, dE))
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

  # Set initial conditions
  initial_conditions1 <- c(Dm2 = 1, E = 0)

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
  LM <- solution2[, "DA1"][[2]]
  logLMb <- log(LM)
  return(logLMb)
}

