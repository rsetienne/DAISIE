#' @name DAISIE_DE_logpNE_unknown_coltime
#' @title Function to calculate the likelihood of observing a non-endemic lineage on the island
#' with unknown colonization time.
#' @description This function calculates the log-likelihood of observing a non-endemic lineage on an island
#' for which the exact colonization time is unknown.
#'
#' @inheritParams default_params_doc_DAISIE_DE
#' @return The output is a numeric value representing the log-likelihood of observing an endemic singleton lineage
#' with unknown colonization time.
#' \item{logL0b}{ The log-likelihood value computed based on system of differential equations.}
#'#' @examples
#'
#' # Select a dataset from a DAISIE package
#'
#' data(Galapagos_datalist)
#' datalist <- Galapagos_datalist
#'
#' # Select a non-endemic lineage in the dataset
#' i <- 2
#' # Define example parameters
#' pars1 <- c(0.2, 0.1, 0.05, 0.02, 0.03)
#'
#' # choose the method to solve the system of differential equations
#' log_likelihood <- DAISIE_DE_logpNE(datalist, i, pars1, methode = "lsodes", reltolint = 1e-16, abstolint = 1e-16)
#'
#' print(log_likelihood)
#'
#' @export DAISIE_DE_logpNE_unknown_coltime


DAISIE_DE_logpNE_unknown_coltime <- function(datalist,
                                             i,
                                             pars1,
                                             methode,
                                             rtol,
                                             atol) {
  t0 <- datalist[[1]]$island_age
  tp <- 0
  parameters <- pars1

  # Define system of equations for interval [t0, tp]
  interval1 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dDA1 <- -pars1[4] * DA1 + pars1[4] * Dm1
      dDm1 <- -(pars1[5] + pars1[1] + pars1[3]) * Dm1 + (pars1[5] * E + pars1[1] * E^2 + pars1[3]) * DA1
      dE1 <- pars1[2] - (pars1[1] + pars1[2]) * E + pars1[1] * E^2
      list(c(dDA1, dDm1, dE))
    })
  }

  # Set initial conditions
  initial_conditions1 <- c(DA1 = 0, Dm1 = 1, E = 0)

  # Time sequence for interval [t0, tp]
  time1 <- c(tp, t0)

  # Solve the system for interval [t0, tp]
  solution1 <- deSOlve::ode(y = initial_conditions1,
                            times = time1,
                            func = interval1,
                            parms = parameters,
                            method = methode,
                            rtol = rtol,
                            atol = atol)

  # Extract log-likelihood
  L0 <- solution1[, "DA1"][[2]]
  logL0b <- log(L0)
  return(logL0b)
}
