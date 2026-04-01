#' @name DAISIE_DE_logpEC_general
#' @title Function to calculate the likelihood of observing an endemic lineage
#' with fixed colonization time. This is valid for infinite K according to the
#' DE equations.
#' @description This function calculates the log-likelihood of observing an
#' endemic lineage with fixed colonization time. This is valid for infinite K
#' according to the DE equations.
#' @inheritParams default_params_doc
#' @return the loglikelihood
#' @examples
#'
#' # Select a dataset from a DAISIE package
#'
#' data(Galapagos_datalist)
#' datalist <- Galapagos_datalist
#' brts <- datalist[[4]]$branching_times
#' missnumspec <- datalist[[4]]$missing_species
#'
#' # Define example parameters
#' pars1 <- c(0.2, 0.1, 0.05, 0.02, 0.03)
#'
#' # choose the method to solve the system of differential equations
#' log_likelihood <- DAISIE_DE_logpEC_general(brts = brts,
#'                                    missnumspec = missnumspec,
#'                                    pars1 = pars1,
#'                                    methode = "lsodes",
#'                                    reltolint = 1e-16,
#'                                    abstolint = 1e-16)
#' @noRd

DAISIE_DE_logpEC_general <- function(brts,
                                     missnumspec,
                                     stac = 0,
                                     pars1,
                                     methode,
                                     reltolint,
                                     abstolint,
                                     use_rcpp = FALSE) {

  if (!(stac %in% c(2, 3, 6))) {
    stop("stac must be 2, 3, or 6 for this function.")
  }


  t0 <- brts[1]
  t1 <- brts[2]
  t2 <- brts[3]
  tp <- 0
  ti <- sort(brts)
  ti <- ti[1:(length(ti) - 2)]

  # Initial conditions
  number_of_species <- length(brts) - 1
  rho <- number_of_species / (missnumspec + number_of_species)

  initial_conditions1   <- c(DE = rho, DM3 = 0, E = 1 - rho, DA3 = 1)
  if (stac == 3) {
    initial_conditions1 <- c(DE = rho, DM3 = 1, E = 1 - rho, DA3 = 0)
  }

  solution0 <- DAISIE_DE_solve_branch(interval_func = interval2_EC,
                                      initial_conditions = initial_conditions1,
                                      time = c(0, ti),
                                      parameter = pars1,
                                      methode = methode,
                                      rtol = reltolint,
                                      atol = abstolint,
                                      use_rcpp = use_rcpp)

  # Time sequences for interval [t2, tp]
  times <- rbind(c(0, ti[1:(length(ti) - 1)]), ti)

  for (idx in 1:length(ti)) {
    # Time sequence idx in interval [t2, tp]
    time1 <- times[, idx]

    # Solve the system for interval [t2, tp]
    solution1 <- DAISIE_DE_solve_branch(interval_func = interval2_EC,
                                        initial_conditions = initial_conditions1,
                                        time = time1,
                                        parameter = pars1,
                                        methode = methode,
                                        rtol = reltolint,
                                        atol = abstolint,
                                        use_rcpp = use_rcpp)

    initial_conditions1 <- c(DE = pars1[1] * solution0[, "DE"][idx + 1] * solution1[, "DE"][2],
                             DM3 = 0,
                             E = solution0[, "E"][idx + 1],
                             DA3 = 1)

  }

  # Initial conditions
  if (stac == 6) {
    initial_conditions2 <- c(DE = initial_conditions1["DE"][[1]],
                             DM1 = 0,
                             DM2 = initial_conditions1["DE"][[1]] * solution0[, "DA3"][length(ti) + 1],
                             DM3 = solution0[, "DM3"][length(ti) + 1],
                             E = initial_conditions1["E"][[1]],
                             DA2 = 0,
                             DA3 = solution0[, "DA3"][length(ti) + 1])
    interval_func <- interval3_ES
  } else {
    initial_conditions2 <- c(DE = initial_conditions1["DE"][[1]],
                             DM2 = initial_conditions1["DE"][[1]] * solution0[, "DA3"][length(ti) + 1],
                             DM3 = solution0[, "DM3"][length(ti) + 1],
                             E = initial_conditions1["E"][[1]],
                             DA3 = solution0[, "DA3"][length(ti) + 1])
    interval_func <- interval2_ES
  }

  # Time sequence for interval [t1, t2]
  time2 <- c(t2, t1)

  # Solve the system for interval [t2, tp]
  solution2 <- DAISIE_DE_solve_branch(interval_func = interval_func,
                                      initial_conditions = initial_conditions2,
                                      time = time2,
                                      parameter = pars1,
                                      methode = methode,
                                      rtol = reltolint,
                                      atol = abstolint,
                                      use_rcpp = use_rcpp)

  # Initial conditions
  if (stac == 6) {
    initial_conditions3 <- c(DA1 = solution2[, "DA2"][[2]],
                             DM1 = solution2[, "DM1"][[2]],
                             E   = solution2[, "E"][[2]])
  } else {
    initial_conditions3 <- c(DA1 = pars1[4] * solution2[, "DM2"][[2]],
                             DM1 = pars1[4] * solution2[, "DM2"][[2]],
                             E   = solution2[, "E"][[2]])
  }

  # Time sequence for interval [t0, t1]
  time3 <- c(t1, t0)

  # Solve the system for interval [t0, t1]
  solution3 <- DAISIE_DE_solve_branch(interval_func = interval4,
                                      initial_conditions = initial_conditions3,
                                      time = time3,
                                      parameter = pars1,
                                      methode = methode,
                                      rtol = reltolint,
                                      atol = abstolint,
                                      use_rcpp = use_rcpp)

  # Extract log-likelihood
  Lk <- solution3[, "DA1"][[2]]
  logLkb <- log(Lk)
  return(logLkb)
}
