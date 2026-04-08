#' @name DAISIE_DE_logpES
#' @title Function to calculate the likelihood of observing an endemic singleton
#' lineage with fixed colonization time. This is valid for infinite K according
#' to the DE equations.
#' @description This function calculates the log-likelihood of observing an
#' endemic singleton lineage with fixed colonization time. This is valid for
#' infinite K according to the DE equations.
#' @inheritParams default_params_doc
#' @return the loglikelihood
#' @examples
#'
#' # Select a dataset from a DAISIE package
#'
#' data(Galapagos_datalist)
#' datalist <- Galapagos_datalist
#' brts <- datalist[[9]]$branching_times
#' missnumspec <- datalist[[9]]$missing_species
#' # Define example parameters
#' pars1 <- c(2.546591, 2.678781, 2.678781, 0.009326754, 1.008583)
#'
#' # choose the method to solve the system of differential equations
#' log_likelihood <- DAISIE_DE_logpES(brts = brts,
#'                                   missnumspec = missnumspec,
#'                                    stac = 2,
#'                                    pars1 = pars1,
#'                                    methode = "lsodes",
#'                                    reltolint = 1e-16,
#'                                    abstolint = 1e-16)
#' @noRd
DAISIE_DE_logpES <- function(brts,
                             missnumspec = 0,
                             pars1,
                             stac = 0,
                             methode = "odeint::runge_kutta_cash_karp54",
                             reltolint = 1e-15,
                             abstolint = 1e-15) {


  if (!(stac %in% c(2, 3, 5, 9))) {
    stop("stac must be 2, 3, 5 or 9 for this function.")
  }

  t0 <- brts[1]
  t1 <- brts[2]
  t2 <- brts[3]
  tp <- 0
  parameters <- pars1

  # Initial conditions

  number_of_species <- length(brts) - 1
  rho <- number_of_species / (missnumspec + number_of_species)

  #pES
  initial_conditions1   <- c(DE = rho, DM2 = 0, DM3 = 0, E = 1 - rho, DA3 = 1)
  interval_func = ifelse(startsWith(methode, "odeint::"), "interval2_ES", interval2_ES)
  time1 <- c(tp, t1)
  # mainland
  if (stac == 3) {
    initial_conditions1 <- c(DE = rho, DM2 = 0, DM3 = 1, E = 1 - rho, DA3 = 0)
  } else if (stac == 5) {
    initial_conditions1 <- c(DE = rho, DM1 = 0, DM2 = 0, DM3 = 0, E = 1 - rho, DA2 = 0, DA3 = 1)
    interval_func <- ifelse(startsWith(methode, "odeint::"), "interval3_ES", interval3_ES)
  } else if (stac == 9) {
    initial_conditions1 <- c(DE = 1  , DM2 = 0, DM3 = 0, E = 0, DA3 = 1)
    time1 <- c(tp, t2)
  }


  # Time sequence for interval [t1, tp]
  solution1 <- DAISIE_DE_solve_branch(interval_func = interval_func,
                                      initial_conditions = initial_conditions1,
                                      parameter = parameters,
                                      time = time1,
                                      methode = methode,
                                      atol = abstolint,
                                      rtol = reltolint)


  if (stac == 9) {
    initial_conditions2 <- c(DE  = solution1[, "DE"][[2]],
                             DM1 = 0,
                             DM2 = solution1[, "DM2"][[2]],
                             DM3 =  solution1[, "DM3"][[2]],
                             E   =  solution1[, "E"][[2]],
                             DA2 = 0,
                             DA3 = solution1[, "DA3"][[2]])
    # Time sequence for interval [t2, t1]
    time2 <- c(t2, t1)

    # Solve the system for interval [t2, t1]
    solution2 <- DAISIE_DE_solve_branch(interval_func = interval3_ES,
                                        initial_conditions = initial_conditions2,
                                        time = time2,
                                        parameter = parameters,
                                        methode = methode,
                                        rtol = reltolint,
                                        atol = abstolint)

  }

  # Initial conditions
  initial_conditions3 <- c(DA1 = pars1[4] * solution1[, "DM2"][[2]],
                           DM1 = pars1[4] * solution1[, "DM2"][[2]],
                           E   = solution1[, "E"][[2]])

  if (stac == 9) {
    initial_conditions3 <- c(DA1 = solution2[, "DA2"][[2]],
                             DM1 = solution2[, "DM1"][[2]],
                             E   = solution2[, "E"][[2]])
  } else if (stac == 5) {
    initial_conditions3 <- c(DA1 = solution1[, "DA2"][[2]],
                             DM1 = solution1[, "DM1"][[2]],
                             E   = solution1[, "E"][[2]])
  }

  # Time sequence for interval [t0, t1]
  time2 <- c(t1, t0)

  solution2 <- DAISIE_DE_solve_branch(interval_func = interval4,
                                      initial_conditions = initial_conditions3,
                                      parameter = parameters,
                                      time = time2,
                                      methode = methode,
                                      atol = abstolint,
                                      rtol = reltolint)

  # Extract log-likelihood
  L1 <- solution2[, "DA1"][[2]]
  logL1b <- log(L1)
  return(logL1b)
}
