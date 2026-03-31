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
#'                                    missnumspec = missnumspec,
#'                                    pars1 = pars1,
#'                                    methode = "lsodes",
#'                                    reltolint = 1e-15,
#'                                    abstolint = 1e-15)
#' @noRd
DAISIE_DE_logpES <- function(brts,
                             missnumspec,
                             pars1,
                             methode,
                             reltolint,
                             abstolint,
                             use_rcpp = FALSE) {
  t0 <- brts[1]
  t1 <- brts[2]
  tp <- 0
  parameters <- pars1

  # Initial conditions

  number_of_species <- length(brts) - 1
  rho <- number_of_species / (missnumspec + number_of_species)

  initial_conditions1 <- c(DE = rho,  DM2 = 0, DM3 = 0, E = 1 - rho, DA3 = 1)

  # Time sequence for interval [t1, tp]
  time1 <- c(tp, t1)

  solution1 <- DAISIE_DE_solve_branch(interval_func = interval2_ES,
                                      initial_conditions = initial_conditions1,
                                      parameter = parameters,
                                      time = time1,
                                      methode = methode,
                                      atol = abstolint,
                                      rtol = reltolint,
                                      use_rcpp = use_rcpp)

  # Initial conditions
  initial_conditions2 <- c(DA1 = pars1[4] * solution1[, "DM2"][[2]],
                           DM1 = pars1[4] * solution1[, "DM2"][[2]],
                           E = solution1[, "E"][[2]])

  # Time sequence for interval [t0, t1]
  time2 <- c(t1, t0)

  solution2 <- DAISIE_DE_solve_branch(interval_func = interval4,
                                      initial_conditions = initial_conditions2,
                                      parameter = parameters,
                                      time = time2,
                                      methode = methode,
                                      atol = abstolint,
                                      rtol = reltolint,
                                      use_rcpp = use_rcpp)

  # Extract log-likelihood
  L1 <- solution2[, "DA1"][[2]]
  logL1b <- log(L1)
  return(logL1b)
}



