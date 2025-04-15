#' @name DAISIE_DE_logpES_unknown_coltime
#' @title Function to calculate the likelihood of observing an endemic singleton lineage on the island
#' with unknown colonization time.
#' @description This function calculates the log-likelihood of observing an endemic singleton lineage on an island
#' for which the exact colonization time is unknown.
#'
#' @inheritParams default_params_doc_DAISIE_DE
#' @return The output is a numeric value representing the log-likelihood of observing an endemic singleton lineage
#' with unknown colonization time.
#' \item{logL1b}{ The log-likelihood value computed based on a system of differential equations.}
#'#'
#' @examples
#'
#' # Select a dataset from a DAISIE package
#'
#' data(Biwa_datalist)
#' datalist <- Biwa_datalist
#'
#' # Select an endemic lineage in the dataset
#' i <- 60
#' # Define example parameters
#' pars1 <- c(0.2, 0.1, 0.05, 0.02, 0.03)
#'
#' # choose the method to solve the system of differential equations
#' log_likelihood <- DAISIE_DE_logpES_unknown_coltime(datalist, i, pars1, methode = "lsodes", reltolint = 1e-16, abstolint = 1e-16)
#'
#' print(log_likelihood)


#' @export DAISIE_DE_logpES_unknown_coltime



### Using D-E approach
DAISIE_DE_logpES_unknown_coltime <- function(datalist,
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

      dDE <- -(pars1[1] + pars1[2]) * DE + 2 * pars1[1] * DE * E

      dDA3 <- -pars1[4] * DA3 + pars1[4] * Dm3

      dDA2 <- -pars1[4] * DA2 + pars1[4] * Dm2

      dDm3 <- -(pars1[5] + pars1[1] + pars1[3]) * Dm3 + (pars1[5] * E + pars1[1] * E^2 + pars1[3]) * DA3

      dDm2 <- -(pars1[5] + pars1[1] + pars1[3]) * Dm2 + (pars1[5] * DE + 2 * pars1[1] * DE * E) * DA3 +
        (pars1[3] + pars1[5] * E + pars1[1] * E^2)* DA2

      dE <- pars1[2] - (pars1[1] + pars1[2]) * E + pars1[1] * E^2

      list(c(dDE, dDA3, dDA2, dDm3, dDm2, dE))
    })
  }

  # Initial conditions
  number_of_species <- length(brts) -1
  number_of_missing_species <- missnumspec
  ro <- number_of_species / (number_of_missing_species + number_of_species)

  if (missnumspec== 0)

  {
    initial_conditions1 <- c(DE = 1, DA3 = 1, D02 = 0, Dm3 = 0, Dm2 = 0, E = 0)
  }
  else

  {
    initial_conditions1 <- c(DE = ro, DA3 = 1, DA2 = 0, Dm3 = 0, Dm2 = 0, E = 1 - ro)
  }
  # Time sequence for interval [t1, tp]
  time1 <- c(tp, t0)

  # Solve the system for interval [t1, tp]
  solution1 <- deSolve::ode(y = initial_conditions1,
                            times = time1,
                            func = interval1,
                            parms = parameters,
                            method = methode,
                            rtol = reltolint,
                            atol = abstolint)

  # Extract log-likelihood
  L1 <- solution1[, "DA2"][[2]]
  logL1b <- log(L1)
  return(logL1b)
}
