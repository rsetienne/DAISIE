#' @name DAISIE_DE_logpEC_max_age_coltime_and_mainland
#' @title Function to calculate the likelihood of observing an endemic lineage on the island
#' with maximum age of colonization, and that coexists with its mainland ancestors. This is valid
#' for infinite K according to the DE equations.
#' @description This function calculates the log-likelihood of observing an endemic lineage on an island
#' with maximum age of colonization, and that coexists with its mainland ancestors. This is valid
#' for infinite K according to the DE equations.
#' @inheritParams default_params_doc
#' @return the loglikelihood
#' @noRd

DAISIE_DE_logpEC_max_age_coltime_and_mainland <- function(brts,
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
  ti <- sort(brts)
  ti <- ti[1:(length(ti)-2)]
  parameters <- pars1

  # Initial conditions
  number_of_species <- length(brts) - 1
  rho <- number_of_species / (missnumspec + number_of_species)

  initial_conditions1 <- c(DE = rho, DA3 = 0, Dm3 = 1, E = 1 - rho)

  solution0 <- DAISIE_DE_solve_branch(interval_func = interval1_3,
                             initial_conditions = initial_conditions1,
                             parameter = pars1,
                             time = c(0, ti),
                             methode = methode,
                             atol = abstolint,
                             rtol = reltolint,
                             use_rcpp = use_rcpp)

  # Time sequences for interval [t2, tp]
  times <- rbind(c(0, ti[1:(length(ti) - 1)]), ti)

  for (idx in 1:length(ti)) {
    # Time sequence idx in interval [t2, tp]
    time1 <- times[, idx]

    # Solve the system for interval [t2, tp]
    solution1 <- DAISIE_DE_solve_branch(interval_func = interval1_3,
                             initial_conditions = initial_conditions1,
                             parameter = pars1,
                             time = time1,
                             methode = methode,
                             atol = abstolint,
                             rtol = reltolint,
                             use_rcpp = use_rcpp)

    # Initial conditions
    initial_conditions1 <- c(DE = pars1[1] * solution0[, "DE"][idx + 1] * solution1[, "DE"][2],
                             DA3 = 1, Dm3 = 0, E = solution0[, "E"][idx + 1])
  }

  # Initial conditions
  initial_conditions2 <- c(DE = initial_conditions1["DE"][[1]],
                           DA2 = 0,
                           DA3 = solution0[, "DA3"][length(ti) + 1],
                           Dm1 = 0,
                           Dm2 = initial_conditions1["DE"][[1]] * solution0[, "DA3"][length(ti)+1],
                           Dm3 = solution0[, "Dm3"][length(ti) + 1],
                           E = initial_conditions1["E"][[1]])

  # Time sequence for interval [t1, t2]
  time2 <- c(t2, t1)

  # Solve the system for interval [t2, tp]
  solution2 <- DAISIE_DE_solve_branch(interval_func = interval3,
                             initial_conditions = initial_conditions2,
                             parameter = pars1,
                             time = time2,
                             methode = methode,
                             atol = abstolint,
                             rtol = reltolint,
                             use_rcpp = use_rcpp)

  # Time sequence for interval [t1, tp]
  time3 <- c(t1, t0)

  # Initial conditions
  initial_conditions3 <- c(DA1 = solution2[, "DA2"][[2]],
                           Dm1 = solution2[, "Dm1"][[2]],
                           E = solution2[, "E"][[2]])

  # Solve the system for interval [t0, t1]
  solution3 <- DAISIE_DE_solve_branch(interval_func = interval3,
                             initial_conditions = initial_conditions3,
                             parameter = pars1,
                             time = time3,
                             methode = methode,
                             atol = abstolint,
                             rtol = reltolint,
                             use_rcpp = use_rcpp)

  # Extract log-likelihood
  Lk <- (solution3[, "DA1"][[2]])
  logLkb <- log(Lk)
  return(logLkb)
}
