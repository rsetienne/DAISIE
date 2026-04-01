#' @name DAISIE_DE_logpNE_general
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
#' pars1 <- c(2.546591, 2.678781, 2.678781, 0.009326754, 1.008583)
#'
#' # choose the method to solve the system of differential equations
#' log_likelihood <- DAISIE_DE_logpNE(brts = brts,
#'                                    pars1 = pars1,
#'                                    methode = "lsodes",
#'                                    reltolint = 1e-16,
#'                                    abstolint = 1e-16)
#' @noRd
DAISIE_DE_logpNE_general <- function(brts,
                                     pars1,
                                     stac = 0,
                                     methode,
                                     reltolint,
                                     abstolint,
                                     use_rcpp = FALSE) {

  if (!(stac %in% c(1, 4, 8))) {
    stop("NE only supports stac values of 1, 4 and 8")
  }

  t0 <- brts[1]
  t1 <- brts[2]
  t2 <- brts[3]
  tp <- 0
  parameters <- pars1

  # Set initial conditions
  interval_func = interval2_NE
  initial_conditions1 <- c(DM2 = 1, E = 0)
  if (stac == 1) { # NE_max_age
    interval_func = interval3_NE
    initial_conditions1 <- c(DM1 = 0, DM2 = 1, E = 0, DA2 = 0)
  }

  # Time sequence for interval [t1, tp]
  time1 <- c(tp, t1)

  # Solve the system for interval [t1, tp]
  solution1 <- DAISIE_DE_solve_branch(interval_func = interval_func,
                                      initial_conditions = initial_conditions1,
                                      time = time1,
                                      parameter = parameters,
                                      methode = methode,
                                      rtol = reltolint,
                                      atol = abstolint,
                                      use_rcpp = use_rcpp)


  if (stac == 8) { # max_min age
    initial_conditions2 <- c(DM1 = 0,
                             DM2 = solution1[, "DM2"][[2]],
                             E   = solution1[, "E"][[2]],
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
  }

  time2 <- c(t1, t0)

  if (stac == 1) { #NE_max_age
    initial_conditions2 <- c(DA1 = solution1[, "DA2"][[2]],
                             DM1 = solution1[, "DM1"][[2]],
                             E   = solution1[, "E"][[2]])
  } else if (stac == 8) { #NE_max_min_age
    initial_conditions2 <- c(DA1 = solution2[, "DA2"][[2]],
                             DM1 = solution2[, "DM1"][[2]],
                             E   = solution2[, "E"][[2]])
  } else { # stac = 4, NE
    initial_conditions2 <- c(DA1 = pars1[4] * solution1[, "DM2"][[2]],
                             DM1 = pars1[4] * solution1[, "DM2"][[2]],
                             E   = solution1[, "E"][[2]])
  }

  # Time sequence for interval [t0, t1]


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
  LM <- solution2[, "DA1"][[2]]
  logLMb <- log(LM)
  return(logLMb)
}
