#' @name DAISIE_DE_logp0
#' @title Log-likelihood of a lineage that colonizes the island but leaves no descendants
#'
#' @description
#' Computes the log-likelihood of a colonization event in which a lineage arrives on the island
#' but does not leave any surviving descendants (i.e., it goes extinct without speciation or persistence).
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
#'
#' @param pars1 A numeric vector of model parameters:
#' \itemize{
#'   \item \code{pars1[1]}: \eqn{\lambda^c} (Cladogenesis rate)
#'   \item \code{pars1[2]}: \eqn{\mu_E} (Extinction rate of endemic lineages)
#'   \item \code{pars1[3]}: \eqn{\mu_{NE}} (Extinction rate of non-endemic lineages)
#'   \item \code{pars1[4]}: \eqn{\gamma} (Colonization rate)
#'   \item \code{pars1[5]}: \eqn{\lambda^a} (Anagenesis rate)
#' }
#'
#' @param methode A character string specifying the numerical method to use for solving the ODEs.
#'
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
      dD0 <- -pars1[4] * D0 + pars1[4] * Dm
      dDm <- -(pars1[5] + pars1[1] + pars1[3]) * Dm + (pars1[5] * E1 + pars1[1] * E1^2 + pars1[3]) * D0
      dE1 <- pars1[2] - (pars1[1] + pars1[2]) * E1 + pars1[1] * E1^2
      list(c(dD0, dDm, dE1))
    })
  }

  # Set initial conditions
  initial_conditions0 <- c(D0 = 1, Dm = 0, E1 = 0)

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
  L0 <- solution0[, "D0"][[2]]
  logL0b <- log(L0)
  return(logL0b)
}




