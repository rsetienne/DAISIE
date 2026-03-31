#' @name DAISIE_DE_logpES_max_min_age_coltime
#' @title Function to compute the likelihood of observing an endemic singleton lineage
#' on the island given the minimum and maximum colonization ages, valid for infinite K
#'according to the DE equations.
#' @description This function calculates the log-likelihood of observing an endemic singleton lineage on an island
#' for which the exact colonization time is unknown, but the maximum and minimum ages of colonization are given.
#' This is valid for infinite K according to the DE equations.
#' @inheritParams default_params_doc
#' @return the loglikelihood
#' @noRd

DAISIE_DE_logpES_max_min_age_coltime <- function(brts,
                                                 missnumspec,
                                                 pars1,
                                                 methode,
                                                 reltolint,
                                                 abstolint,
                                                 use_rcpp = FALSE) {
  t0 <- brts[1]
  t1 <- brts[2]
  t2 <- brts[3]
  tp <- 0
  parameters <- pars1

  # Initial conditions

  initial_conditions1 <- c(DE = 1, DM2 = 0, DM3 = 0, E = 0, DA3 = 1)


  # Time sequence for interval [tp, t2]
  time1 <- c(tp, t2)

  # Solve the system for interval [tp, t2]
  solution1 <- DAISIE_DE_solve_branch(interval_func = interval2_ES,
                                      initial_conditions = initial_conditions1,
                                      time = time1,
                                      parameter = parameters,
                                      methode = methode,
                                      rtol = reltolint,
                                      atol = abstolint,
                                      use_rcpp = use_rcpp)
  # Initial conditions

  initial_conditions2 <- c(DE = solution1[, "DE"][[2]],
                           DM1 = 0,
                           DM2 = solution1[, "DM2"][[2]],
                           DM3 =  solution1[, "DM3"][[2]],
                           E =  solution1[, "E"][[2]],
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
                                      atol = abstolint,
                                      use_rcpp = use_rcpp)

  # Initial conditions
  initial_conditions3 <- c(DA1 = solution2[, "DA2"][[2]],
                           DM1 = solution2[, "DM1"][[2]],
                           E = solution2[, "E"][[2]])

  # Time sequence for interval [t1, t0]
  time3 <- c(t1, t0)

  # Solve the system for interval [t1, t0]
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
