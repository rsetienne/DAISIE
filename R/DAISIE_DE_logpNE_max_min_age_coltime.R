#' @name DAISIE_DE_logpNE_max_min_age_coltime
#' @title Function to calculate the likelihood of observing a non-endemic lineage on the island
#' with minimum and maximum times of colonization. This valid for infinite K according to the DE equations.
#' @description This function calculates the log-likelihood of observing a non-endemic lineage on an island
#' for which the exact colonization time is unknown, but the maximum and minimum times of colonization are
#' known. This is valid for infinite K according to the DE equations
#' @inheritParams default_params_doc
#' @return the loglikelihood
#' @noRd

DAISIE_DE_logpNE_max_min_age_coltime <- function(brts,
                                                 pars1,
                                                 methode = "ode45",
                                                 reltolint = 1e-12,
                                                 abstolint = 1e-12,
                                                 use_rcpp = FALSE) {
  t0 <- brts[1]
  t1 <- brts[2]
  t2 <- brts[3]
  tp <- 0
  parameters <- pars1

  # Time sequence for interval [t2, tp]
  time1 <- c(tp, t2)

  # Initial conditions
  initial_conditions1 <- c(DM2 = 1, E = 0)

  # Solve the system for interval [t2, tp]
  solution1 <- DAISIE_DE_solve_branch(interval_func = interval2_NE,
                                      initial_conditions = initial_conditions1,
                                      time = time1,
                                      parameter = parameters,
                                      methode = methode,
                                      rtol = reltolint,
                                      atol = abstolint,
                                      use_rcpp = use_rcpp)

  # Initial conditions
  initial_conditions2 <- c(DM1 = 0,
                           DM2 = solution1[, "DM2"][[2]],
                           E = solution1[, "E"][[2]],
                           DA2 = 0)

  # Time sequence for interval [t1, t2]
  time2 <- c(t2, t1)

  # Solve the system for interval [t1, tp]
  solution2 <- DAISIE_DE_solve_branch(interval_func = interval3_NE,
                                      initial_conditions = initial_conditions2,
                                      time = time2,
                                      parameter = parameters,
                                      methode = methode,
                                      rtol = reltolint,
                                      atol = abstolint,
                                      use_rcpp = use_rcpp)

  # Initial conditions
  initial_conditions3 <- c(DA1 = solution2[, "DA2"][[2]],
                           DM1 = solution2[, "DM1"][[2]],
                           E = solution2[, "E"][[2]])

  # Time sequence for interval [t0, t1]
  time3 <- c(t1, t0)

  # Solve the system for interval [t0, t1]
  solution3 <- DAISIE_DE_solve_branch(interval_func = interval4,
                                      initial_conditions = initial_conditions3,
                                      time = time3,
                                      parameter = parameters,
                                      methode = methode,
                                      rtol = reltolint,
                                      atol = abstolint,
                                      use_rcpp = use_rcpp)

  # Extract log-likelihood
  L1 <- solution3[, "DA1"][[2]]
  logL1b <- log(L1)
  return(logL1b)
}
