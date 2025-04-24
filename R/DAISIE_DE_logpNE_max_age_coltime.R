#' @name DAISIE_DE_logpNE_max_age_coltime
#' @title Function to calculate the likelihood of observing a non-endemic lineage on the island
#' with a maximum age of colonization.
#' @description This function calculates the log-likelihood of observing a non-endemic lineage on an island
#' for which the exact colonization time is unknown, but the maximum age of colonization is given.
#'
#' @inheritParams default_params_doc_DAISIE_DE
#'
#' @return The output is a numeric value representing the log-likelihood of observing a non-endemic singleton lineage
#' for which the maximum colonization time is given.
#' \item{logL1b}{ The log-likelihood value computed based on a system of differential equations. }
#' @examples
#'
#' # Select a dataset from a DAISIE package
#'
#' data(Biwa_datalist)
#' datalist <- Biwa_datalist
#'
#' # Select an endemic lineage in the dataset
#' i <- 40
#' # Define example parameters
#' pars1 <- c(0.2, 0.1, 0.05, 0.02, 0.03)
#'
#' # choose the method to solve the system of differential equations
#' log_likelihood <- DAISIE_DE_logpNE_max_age_coltime(datalist, i, pars1, methode = "lsodes", reltolint = 1e-16, abstolint = 1e-16)
#'
#' print(log_likelihood)
#'


#' @export DAISIE_DE_logpNE_max_age_coltime



### Using D-E approach
DAISIE_DE_logpNE_max_age_coltime <- function(datalist,
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

  interval1 <- function(t, state, parameters) {
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
  interval2 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dDA <- -pars1[4] * DA + pars1[4] * Dm1
      dDm1 <- -(pars1[5] + pars1[1] + pars1[3]) * Dm1 + (pars1[5] * E + pars1[1] * E^2 + pars1[3]) * DA
      dE <- pars1[2] - (pars1[1] + pars1[2]) * E + pars1[1] * E^2
      list(c(dDA, dDm1, dE))
    })
  }

  # Initial conditions
  initial_conditions1 <- c(DA = 0, Dm1 = 0, Dm2 = 1, E = 0)

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
  initial_conditions2 <- c(DA = solution1[, "DA"][[2]],
                           Dm1 = solution1[, "Dm1"][[2]],
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
  L1 <- solution2[, "DA"][[2]]
  logL1b <- log(L1)
  return(logL1b)
}



