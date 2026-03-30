#' @name DAISIE_DE_logpNE_max_age_coltime
#' @title Function to calculate the likelihood of observing a non-endemic
#' lineage on the island with a maximum time of colonization. This is valid for
#' infinite K according to the DE equations.
#' @description This function calculates the log-likelihood of observing a
#' non-endemic lineage on an island for which the exact colonization time is
#' unknown, but the maximum time of colonization is given. This is valid for
#' infinite K according to the DE equations.
#' @inheritParams default_params_doc
#' @return the loglikelihood
#' @examples
#'
#' # Select a dataset from a DAISIE package
#'
#' data(Biwa_datalist)
#' datalist <- Biwa_datalist
#' brts <- datalist[[40]]$branching_times
#' # Define example parameters
#' pars1 <- c(0.2, 0.1, 0.05, 0.02, 0.03)
#'
#' # choose the method to solve the system of differential equations
#' log_likelihood <- DAISIE_DE_logpNE_max_age_coltime(brts = brts,
#'                                                    pars1 = pars1,
#'                                                    methode = "lsodes",
#'                                                    reltolint = 1e-16,
#'                                                    abstolint = 1e-16)
#' @noRd

DAISIE_DE_logpNE_max_age_coltime <- function(brts,
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
  initial_conditions1 <- c(DA = 0, Dm1 = 0, Dm2 = 1, E = 0)

  # Time sequence for interval [t1, tp]
  time1 <- c(tp, t1)

  solution1 <- DAISIE_DE_solve_branch(interval_func = interval1,
                                      initial_conditions = initial_conditions1,
                                      time = time1,
                                      parameter = parameters,
                                      methode = methode,
                                      rtol = reltolint,
                                      atol = abstolint,
                                      use_rcpp = use_rcpp)

  # Initial conditions
  initial_conditions2 <- c(DA = solution1[, "DA"][[2]],
                           Dm1 = solution1[, "Dm1"][[2]],
                           E = solution1[, "E"][[2]])

  # Time sequence for interval [t0, t1]
  time2 <- c(t1, t0)

  # Solve the system for interval [t0, t1]
  solution2 <- DAISIE_DE_solve_branch(interval_func = interval3,
                                      initial_conditions = initial_conditions2,
                                      time = time2,
                                      parameter = parameters,
                                      methode = methode,
                                      rtol = reltolint,
                                      atol = abstolint,
                                      use_rcpp = use_rcpp)

  # Extract log-likelihood
  L1 <- solution2[, "DA"][[2]]
  logL1b <- log(L1)
  return(logL1b)
}
