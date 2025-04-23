#' @name DAISIE_DE_logp0
#' @title Log-likelihood of a mainland species that colonizes the island
#' but leaves no descendants
#'
#' @description
#' Computes the log-likelihood of a colonization event where a mainland species
#' arrives on the island but does not leave any surviving descendants.
#'
#' @inheritParams default_params_doc_DAISIE_DE
#'
#' @return A single numeric value:
#' \describe{
#'   \item{logL0b}{Log-likelihood of the scenario where the lineage colonizes the island
#'   but does not leave any surviving descendants, computed from the ODE solution.}
#' }
#'
#' @examples
#' # Example model parameters
#' pars1 <- c(0.2, 0.1, 0.05, 0.02, 0.03)
#'
#' # Compute log-likelihood
#' log_likelihood <- DAISIE_DE_logp0(datalist, pars1, methode = "lsodes")
#' print(log_likelihood)
#'
#' @export


# Define system of equations for interval [t0, tp]


DAISIE_DE_logp0 <- function(datalist,
                            pars1,
                            methode) {
  t0 <- datalist[[1]]$island_age
  tp <- 0
  parameters <- pars1

  interval0 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dDA1 <- -pars1[4] * DA1 + pars1[4] * Dm1
      dDm1 <- -(pars1[5] + pars1[1] + pars1[3]) * Dm1 + (pars1[5] * E + pars1[1] * E^2 + pars1[3]) * DA1
      dE <- pars1[2] - (pars1[1] + pars1[2]) * E + pars1[1] * E^2
      list(c(dDA1, dDm1, dE))
    })
  }

  # Set initial conditions
  initial_conditions0 <- c(DA1 = 1, Dm1 = 0, E = 0)

  # Time sequence for interval [t0, tp]
  time0 <- c(tp, t0)

  # Solve the system for interval [t0, tp]
  solution0 <- deSolve::ode(y = initial_conditions0,
                            times = time0,
                            func = interval0,
                            parms = parameters,
                            method = methode,
                            rtol = 1E-12,
                            atol = 1E-12)

  # Extract log-likelihood
  L0 <- solution0[, "DA1"][[2]]
  logL0b <- log(L0)
  return(logL0b)
}

