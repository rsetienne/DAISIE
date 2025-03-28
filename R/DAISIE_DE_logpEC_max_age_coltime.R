#' @name DAISIE_DE_logpEC_max_age_coltime
#' @title Function to calculate the likelihood of observing an endemic lineage on the island
#' with its mainland ancestor
#' @description This function calculates the log-likelihood of observing an endemic lineage on an island
#' that coexist with its mailand ancestors.
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
#' @return The output is a numeric value representing the log-likelihood of observing an endemic lineage
#' withits mainland ancestors
#' \item{logLkb}{ The log-likelihood value computed based on the differential equation system.}
#'
#' @export DAISIE_DE_logpEC_max_age_coltime



### Using D-E approach
DAISIE_DE_logpEC_max_age_coltime <- function(datalist,
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
  ti <- sort(datalist[[i]]$branching_times)
  ti <- ti[1:(length(ti)-2)]
  parameters <- pars1


  # Define system of equations for interval [t2, tp]
  interval1 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {

      dD1 <- -(pars1[1] + pars1[2]) * D1 + 2 * pars1[1] * D1 * E1

      dD0 <- -pars1[4] * D0 + pars1[4] * Dm3

      dDm3 <- -(pars1[5] + pars1[1] + pars1[3]) * Dm3 + (pars1[5] * E1 + pars1[1] * E1^2 + pars1[3]) * D0

      dE1 <- pars1[2] - (pars1[1] + pars1[2]) * E1 + pars1[1] * E1^2

      list(c(dD1, dD0, dDm3, dE1))
    })
  }

  # Initial conditions
  number_of_species <- length(brts) -1
  number_of_missing_species <- missnumspec
  ro <- number_of_species / (number_of_missing_species + number_of_species)

  initial_conditions1 <- c(D1 = ro, D0 = 1, Dm3 = 0, E1 = 1 - ro)



  # Define system of equations for interval [t1, t2]
  interval2 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {

      dD1 <- -(pars1[1] + pars1[2]) * D1 + 2 * pars1[1] * D1 * E1

      dD02 <- -pars1[4] * D02 + pars1[4] * Dm2

      dD03 <- -pars1[4] * D03 + pars1[4] * Dm3

      dDm1 <- -(pars1[5] + pars1[1] + pars1[3] + pars1[4]) * Dm1 +
        (pars1[3] + pars1[5] * E1 + pars1[1] * E1^2)* D02 + pars1[4] * Dm2

      dDm2 <- -(pars1[5] + pars1[1] + pars1[3]) * Dm2 + (pars1[3] + pars1[5] * E1 + pars1[1] * E1^2)* D02 +
        (pars1[5] * D1 + 2 * pars1[1] * D1 * E1 ) * D03


      dDm3 <- -(pars1[5] + pars1[1] + pars1[3]) * Dm3 + (pars1[3] + pars1[5] * E1 + pars1[1] * E1^2) * D03



      dE1 <- pars1[2] - (pars1[1] + pars1[2]) * E1 + pars1[1] * E1^2


      list(c(dD1, dD02, dD03, dDm1, dDm2, dDm3, dE1))
    })
  }




  interval3 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {

      dD0 <- -pars1[4] * D0 + pars1[4] * Dm1

      dDm1 <- -(pars1[5] + pars1[1] + pars1[3]) * Dm1 + (pars1[5] * E1 + pars1[1] * E1^2 + pars1[3]) * D0

      dE1 <- pars1[2] - (pars1[1] + pars1[2]) * E1 + pars1[1] * E1^2

      list(c(dD0, dDm1, dE1))
    })
  }


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
                             D0 = 1, Dm3 = 0, E1 = solution0[, "E1"][idx + 1])
  }

  # Initial conditions
  initial_conditions2 <- c(D1 = initial_conditions1["D1"][[1]],
                           D02 = 0,
                           D03 = solution0[, "D0"][length(ti) + 1],
                           Dm1 = 0,
                           Dm2 = initial_conditions1["D1"][[1]] * solution0[, "D0"][length(ti)+1],
                           Dm3 = solution0[, "Dm3"][length(ti) + 1],
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




  # Time sequence for interval [t1, tp]
  time3 <- c(t1, t0)

  # Initial conditions
  initial_conditions3 <- c(D0 = solution2[, "D02"][[2]],
                           Dm1 = solution2[, "Dm1"][[2]],
                           E1 = solution2[, "E1"][[2]])

  # Solve the system for interval [t0, t1]
  solution3 <- deSolve::ode(y = initial_conditions3,
                            times = time3,
                            func = interval3,
                            parms = pars1,
                            method = methode,
                            rtol = reltolint,
                            atol = abstolint)



  # Extract log-likelihood
  Lk <- (solution3[, "D0"][[2]])
  logLkb <- log(Lk)
  return(logLkb)
}

