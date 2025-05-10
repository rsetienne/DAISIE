#' @name DAISIE_DE_logpES_max_min_age_coltime
#' @title Function to compute the likelihood of observing an endemic singleton lineage
#' on the island given the minimum and maximum colonization ages.

#' @description This function calculates the log-likelihood of observing an endemic singleton lineage on an island
#' for which the exact colonization time is unknown, but the maximum and minimum ages of colonization are given.
#' @inheritParams default_params_doc
#' @return The output is a numeric value representing the log-likelihood of observing an endemic singleton lineage
#' for which the minimum and maximum ages of colonization are given.
#' \item{logL1b}{ The log-likelihood value computed based on a system of differential equations.}
#'
#' @export DAISIE_DE_logpES_max_min_age_coltime

DAISIE_DE_logpES_max_min_age_coltime <- function(brts,
                                                 missnumspec,
                                                 pars1,
                                                 methode,
                                                 reltolint,
                                                 abstolint) {
  t0 <- brts[1]
  t1 <- brts[2]
  t2 <- brts[3]
  tp <- 0
  parameters <- pars1



  # Define system of equations for interval [tp, t3]
  interval1 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {

      dDE <- -(pars1[1] + pars1[2]) * DE + 2 * pars1[1] * DE * E

      dDA3 <- -pars1[4] * DA3 + pars1[4] * Dm3

      dDm2 <- -(pars1[5] + pars1[1] + pars1[3] + pars1[4]) * Dm2 +
        (pars1[5] * DE + 2 * pars1[1] * DE * E) * DA3



      dDm3 <- -(pars1[5] + pars1[1] + pars1[3]) * Dm3 +
        (pars1[3] + pars1[5] * E + pars1[1] * E^2) * DA3


      dE <- pars1[2] - (pars1[1] + pars1[2]) * E + pars1[1] * E^2

      list(c(dDE, dDA3, dDm2, dDm3, dE))
    })
  }


  # Initial conditions

    initial_conditions1 <- c(DE = 1, DA3 = 1, Dm2 = 0, Dm3 = 0, E = 0)


  # Time sequence for interval [tp, t2]
  time1 <- c(tp, t2)

  # Solve the system for interval [tp, t2]
  solution1 <- deSolve::ode(y = initial_conditions1,
                            times = time1,
                            func = interval1,
                            parms = parameters,
                            method = methode,
                            rtol = reltolint,
                            atol = abstolint)



  # Define system of equations for interval [tp, t2]
  interval2 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {

      dDE <- -(pars1[1] + pars1[2]) * DE + 2 * pars1[1] * DE * E


      dDA2 <- -pars1[4] * DA2 + pars1[4] * Dm2

      dDA3 <- -pars1[4] * DA3 + pars1[4] * Dm3


      dDm1 <- -(pars1[5] + pars1[1] + pars1[3] + pars1[4]) * Dm1 +
        (pars1[3] + pars1[5] * E + pars1[1] * E^2)* DA2 + pars1[4] * Dm2



      dDm2 <- -(pars1[5] + pars1[1] + pars1[3]) * Dm2 +
        (pars1[3] + pars1[5] * E + pars1[1] * E^2)* DA2 +
        (pars1[5] * DE + 2 * pars1[1] * DE * E ) * DA3



      dDm3 <- -(pars1[5] + pars1[1] + pars1[3]) * Dm3 +
        (pars1[3] + pars1[5] * E + pars1[1] * E^2) * DA3


      dE <- pars1[2] - (pars1[1] + pars1[2]) * E + pars1[1] * E^2

      list(c(dDE, dDA2, dDA3, dDm1, dDm2, dDm3, dE))
    })
  }


  # Initial conditions

    initial_conditions2 <- c(DE = solution1[, "DE"][[2]],
                             DA2 = 0,
                             DA3 = solution1[, "DA3"][[2]],
                             Dm1 = 0,
                             Dm2 = solution1[, "Dm2"][[2]],
                             Dm3 =  solution1[, "Dm3"][[2]],
                             E =  solution1[, "E"][[2]])




  # Define system of equations for interval [t0, t1]
  interval3 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {

      dDA1 <- -pars1[4] * DA1 + pars1[4] * Dm1

      dDm1 <- -(pars1[5] + pars1[1] + pars1[3]) * Dm1 + (pars1[5] * E + pars1[1] * E^2 + pars1[3]) * DA1

      dE <- pars1[2] - (pars1[1] + pars1[2]) * E + pars1[1] * E^2

      list(c(dDA1, dDm1, dE))
    })
  }



  # Time sequence for interval [t2, t1]
  time2 <- c(t2, t1)

  # Solve the system for interval [t2, t1]
  solution2 <- deSolve::ode(y = initial_conditions2,
                            times = time2,
                            func = interval2,
                            parms = parameters,
                            method = methode,
                            rtol = reltolint,
                            atol = abstolint)

  # Initial conditions
  initial_conditions3 <- c(DA1 = solution2[, "DA2"][[2]],
                           Dm1 = solution2[, "Dm1"][[2]],
                           E = solution2[, "E"][[2]])

  # Time sequence for interval [t1, t0]
  time3 <- c(t1, t0)

  # Solve the system for interval [t1, t0]
  solution3 <- deSolve::ode(y = initial_conditions3,
                            times = time3,
                            func = interval3,
                            parms = parameters,
                            method = methode,
                            rtol = reltolint,
                            atol = abstolint)

  # Extract log-likelihood
  L1 <- solution3[, "DA1"][[2]]
  logL1b <- log(L1)

  return(logL1b)

}



