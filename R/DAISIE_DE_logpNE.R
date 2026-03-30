#' @name DAISIE_DE_logpNE
#' @title Function to calculate the likelihood of observing a non-endemic lineage
#' with fixed colonization time. This is valid for infinite K according to the DE
#' equations.
#' @description This function calculates the log-likelihood of observing a non-endemic lineage
#' with fixed colonization time. This is valid for infinite K according to the DE
#' equations.
#' @inheritParams default_params_doc
#' @return the loglikelihood
#' @examples
#'
#' # Select a dataset from a DAISIE package
#'
#' data(Galapagos_datalist)
#' datalist <- Galapagos_datalist
#' brts <- datalist[[3]]$branching_times
#' # Define example parameters
#' pars1 <- c(0.2, 0.1, 0.05, 0.02, 0.03)
#'
#' # choose the method to solve the system of differential equations
#' log_likelihood <- DAISIE_DE_logpNE(brts = brts,
#'                                    pars1 = pars1,
#'                                    methode = "lsodes",
#'                                    reltolint = 1e-16,
#'                                    abstolint = 1e-16)
#' @noRd

DAISIE_DE_logpNE <- function(brts,
                             pars1,
                             methode,
                             reltolint,
                             abstolint,
                             use_rcpp = FALSE) {

  t0 <- brts[1]
  t1 <- brts[2]
  tp <- 0
  parameters <- pars1
  # Set initial conditions
  initial_conditions1 <- c(Dm2 = 1, E = 0)

  # Time sequence for interval [t1, tp]
  time1 <- c(tp, t1)

  # Solve the system for interval [t1, tp]
  solution1 <- DAISIE_DE_solve_branch(interval_func = interval1_13,
                                      initial_conditions = initial_conditions1,
                                      time = time1,
                                      parameter = parameters,
                                      methode = methode,
                                      rtol = reltolint,
                                      atol = abstolint,
                                      use_rcpp = use_rcpp)

  # Set initial conditions
  initial_conditions2 <- c(DA1 = pars1[4] * solution1[, "Dm2"][[2]],
                           Dm1 = pars1[4] * solution1[, "Dm2"][[2]],
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
  LM <- solution2[, "DA1"][[2]]
  logLMb <- log(LM)
  return(logLMb)
}

