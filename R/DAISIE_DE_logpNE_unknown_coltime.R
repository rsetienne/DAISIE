#' @name DAISIE_DE_logpNE_unknown_coltime
#' @title Function to calculate the likelihood of observing a non-endemic lineage on the island
#' with unknown colonization time
#' @description This function calculates the log-likelihood of observing a non-endemic lineage on an island
#' for which the exact colonization time is unknowned.
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
#' @return The output is a numeric value representing the log-likelihood of observing an endemic singleton lineage
#' with unknown colonization time.
#' \item{logL0b}{ The log-likelihood value computed based on the differential equation system.}
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
#' log_likelihood <- DAISIE_DE_logpNE(datalist, brts = datalist[[i]]$branching_times, pars1, missnumspec = datalist[[i]]$missing_species, methode = "lsodes", reltolint = 1e-16, abstolint = 1e-16)
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
      dD0 <- -pars1[4] * D0 + pars1[4] * Dm
      dDm <- -(pars1[5] + pars1[1] + pars1[3]) * Dm + (pars1[5] * E1 + pars1[1] * E1^2 + pars1[3]) * D0
      dE1 <- pars1[2] - (pars1[1] + pars1[2]) * E1 + pars1[1] * E1^2
      list(c(dD0, dDm, dE1))
    })
  }

  # Set initial conditions
  initial_conditions1 <- c(D0 = 0, Dm = 1, E1 = 0)

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
  L0 <- solution1[, "D0"][[2]]
  logL0b <- log(L0)
  return(logL0b)
}
