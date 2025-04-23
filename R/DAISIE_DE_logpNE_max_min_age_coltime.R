#' @name DAISIE_DE_logpNE_max_min_age_coltime
#' @title Function to calculate the likelihood of observing a non-endemic lineage on the island
#' with minimun and maximun age of colonization
#' @description This function calculates the log-likelihood of observing a non-endemic lineage on an island
#' for which the exact colonization time is unknown, but the maximum and minimum ages of colonization are known.
#'
#' @param datalist A list containing colonization and branching information for island lineages.
#' This object can be created using the \code{DAISIE_dataprep()} function, or manually constructed.
#' It should be a list with the following structure:
#' \itemize{
#'   \item \code{datalist[[1]]$island_age}: Age of the island.
#'   \item \code{datalist[[1]]$not_present} or (for trait-dependent cases)
#'         \code{datalist[[1]]$not_present_type1} and \code{datalist[[1]]$not_present_type2}:
#'         Number of mainland species not present on the island.
#'   \item Each subsequent element of the list corresponds to a single colonist lineage and includes:
#'     \itemize{
#'       \item \code{$colonist_name}: Name of the species or clade.
#'       \item \code{$branching_times}: A numeric vector starting with the island age, followed by colonization and speciation times.
#'       \item \code{$stac}: Colonist status, one of the following:
#'         \enumerate{
#'           \item Non_endemic_MaxAge: 1
#'           \item Endemic: 2
#'           \item Endemic & Non_Endemic: 3
#'           \item Non_Endemic: 4
#'           \item Endemic_Singleton_MaxAge: 5
#'           \item Endemic_Clade_MaxAge: 6
#'           \item Endemic & Non_Endemic_Clade_MaxAge: 7
#'         }
#'       \item \code{$missing_species}: Number of missing species for the clade (applies to endemic clades only).
#'       \item \code{$type1or2}: Lineage type (1 or 2), used in trait-dependent models.
#'     }
#' }
#' @param brts The branching times of the lineage being considered in the dataset.
#' @param missnumspec The number of missing species in the lineage being considered.
#' @param pars1 A numeric vector of model parameters:
#' \itemize{
#'   \item \code{pars1[1]}: \eqn{\lambda^c} (Cladogenesis rate)
#'   \item \code{pars1[2]}: \eqn{\mu_E} (Extinction rate of endemic lineages)
#'   \item \code{pars1[3]}: \eqn{\mu_{NE}} (Extinction rate of non-endemic lineages)
#'   \item \code{pars1[4]}: \eqn{\gamma} (Colonization rate)
#'   \item \code{pars1[5]}: \eqn{\lambda^a} (Anagenesis rate)
#' }
#'
#' @param methode The numerical method to use for solving the system of differential equations.
#' @param reltolint Relative tolerance for numerical integration.
#' @param abstolint Absolute tolerance for numerical integration.
#'
#' @return The output is a numeric value representing the log-likelihood of observing a non-endemic singleton lineage
#' with minimum and maximum age of colonization
#' \item{logL1b}{ The log-likelihood value computed based on the differential equation system.}
#'
#' @export DAISIE_DE_logpNE_max_min_age_coltime



### Using D-E approach
DAISIE_DE_logpNE_max_min_age_coltime <- function(datalist,
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
  parameters <- pars1

  interval1 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {

      dDM <- -(pars1[5] + pars1[1] + pars1[3] + pars1[4]) * DM

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

      dDm <- -(pars1[5] + pars1[1] + pars1[3] + pars1[4]) * Dm +
        (pars1[5] * E1 + pars1[1] * E1^2 + pars1[3]) * D0 + pars1[4] * DM

      dDM <- -(pars1[5] + pars1[1] + pars1[3]) * DM +
        (pars1[5] * E1 + pars1[1] * E1^2 + pars1[3]) * D0

      dE1 <-  pars1[2] - (pars1[1] + pars1[2]) * E1 + pars1[1] * E1^2

      list(c(dD0, dDm, dDM, dE1))
    })
  }



  # Define system of equations for interval [t0, t1]
  interval3 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dD0 <- -pars1[4] * D0 + pars1[4] * Dm

      dDm <- -(pars1[5] + pars1[1] + pars1[3]) * Dm + (pars1[5] * E1 + pars1[1] * E1^2 + pars1[3]) * D0

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





