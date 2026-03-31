#' @name DAISIE_DE_logpES_mainland
#' @title Function to calculate the likelihood of observing an endemic singleton
#' lineage on the island that coexists with its mainland ancestor. This is valid
#' for infinite K according to the DE equations.
#' @description This function calculates the log-likelihood of observing an
#' endemic singleton lineage on an island that coexists with its mainland
#' ancestor. This is valid for infinite K according to the DE equations.
#' @inheritParams default_params_doc
#' @return the loglikelihood
#' @noRd

DAISIE_DE_logpES_mainland <- function(brts,
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

  number_of_species <- length(brts) - 1
  rho <- number_of_species / (missnumspec + number_of_species)

  # Initial conditions
  initial_conditions1 <- c(DE = rho, DA3 = 0, DM3 = 1, DM2 = 0, E = 1 - rho)

    # Time sequence for interval [t1, tp]
  time1 <- c(tp, t1)

  # Solve the system for interval [t1, tp]
  solution1 <- DAISIE_DE_solve_branch(interval_func = interval2_ES,
                                        initial_conditions = initial_conditions1,
                                        time = time1,
                                        parameter = parameters,
                                        methode = methode,
                                        rtol = reltolint,
                                        atol = abstolint,
                                        use_rcpp = use_rcpp)

  # Initial conditions
  initial_conditions2 <- c(DA1 = pars1[4] * solution1[, "DM2"][[2]],
                           DM1 = pars1[4] * solution1[, "DM2"][[2]],
                           E = solution1[, "E"][[2]])

  # Time sequence for interval [t0, t1]
  time2 <- c(t1, t0)

  # Solve the system for interval [t0, t1]
  solution2 <- DAISIE_DE_solve_branch(interval_func = interval4,
                                        initial_conditions = initial_conditions2,
                                        time = time2,
                                        parameter = parameters,
                                        methode = methode,
                                        rtol = reltolint,
                                        atol = abstolint,
                                        use_rcpp = use_rcpp)

  # Extract log-likelihood
  L1 <- solution2[, "DA1"][[2]]
  logL1b <- log(L1)
  return(logL1b)
}
