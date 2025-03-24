###############################################################################
### function to calculate the likelihood of observing an endemic singleton lineage
### with the max and min age colonization time t1
###############################################################################
### Using D-E approach


# pars1[1] corresponds to the Cladogenesis rate
# pars1[2] corresponds to the Extinction rate of endemic lineages
# pars1[3] corresponds to the Extinction rate of non-endemic lineages
# pars1[4] = corresponds to the Colonization rate
# pars1[5] = corresponds to the Anagenesis rate

### Using D-E approach
DAISIE_DE_logpES_max_min_age_coltime <- function(datalist,
                                             i,
                                             pars1,
                                             methode,
                                             reltolint,
                                             abstolint) {

  t0 <- datalist[[i]]$branching_times[1]
  t1 <- datalist[[i]]$branching_times[2]
  t2 <- datalist[[i]]$branching_times[3]
  tp <- 0
  parameters <- pars1



  # Define system of equations for interval [tp, t3]
  interval1 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {

      dD1 <- -(pars1[1] + pars1[2]) * D1 + 2 * pars1[1] * D1 * E1

      dD03 <- -pars1[4] * D03 + pars1[4] * Dm3

      dDm2 <- -(pars1[5] + pars1[1] + pars1[2] + pars1[4]) * Dm2 +
        (pars1[5] * D1 + 2 * pars1[1] * D1 * E1) * D03



      dDm3 <- -(pars1[5] + pars1[1] + pars1[2]) * Dm3 +
        (pars1[2] + pars1[5] * E1 + pars1[1] * E1^2) * D03


      dE1 <- pars1[2] - (pars1[1] + pars1[2]) * E1 + pars1[1] * E1^2

      list(c(dD1, dD03, dDm2, dDm3, dE1))
    })
  }


  # Initial conditions

    initial_conditions1 <- c(D1 = 1, D03 = 1, Dm2 = 0, Dm3 = 0, E1 = 0)


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
      dD1 <- -(pars1[1] + pars1[2]) * D1 + 2 * pars1[1] * D1 * E1


      dD02 <- -pars1[4] * D02 + pars1[4] * Dm2

      dD03 <- -pars1[4] * D03 + pars1[4] * Dm3


      dDm1 <- -(pars1[5] + pars1[1] + pars1[2] + pars1[4]) * Dm1 +
        (pars1[2] + pars1[5] * E1 + pars1[1] * E1^2)* D02 + pars1[4] * Dm2



      dDm2 <- -(pars1[5] + pars1[1] + pars1[2]) * Dm2 +
        (pars1[2] + pars1[5] * E1 + pars1[1] * E1^2)* D02 +
        (pars1[5] * D1 + 2 * pars1[1] * D1 * E1 ) * D03



      dDm3 <- -(pars1[5] + pars1[1] + pars1[2]) * Dm3 +
        (pars1[2] + pars1[5] * E1 + pars1[1] * E1^2) * D03


      dE1 <- pars1[2] - (pars1[1] + pars1[2]) * E1 + pars1[1] * E1^2

      list(c(dD1, dD02, dD03, dDm1, dDm2, dDm3, dE1))
    })
  }


  # Initial conditions

    initial_conditions2 <- c(D1 = solution1[, "D1"][[2]],
                             D02 = 0,
                             D03 = solution1[, "D03"][[2]],
                             Dm1 = 0,
                             Dm2 = solution1[, "Dm2"][[2]]*solution1[, "D03"][[2]],
                             Dm3 =  solution1[, "Dm3"][[2]],
                             E1 =  solution1[, "E1"][[2]])




  # Define system of equations for interval [t0, t1]
  interval3 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {

      dD0 <- -pars1[4] * D0 + pars1[4] * Dm1

      dDm1 <- -(pars1[5] + pars1[1] + pars1[2]) * Dm1 + (pars1[5] * E1 + pars1[1] * E1^2 + pars1[2]) * D0

      dE1 <- pars1[2] - (pars1[1] + pars1[2]) * E1 + pars1[1] * E1^2

      list(c(dD0, dDm1, dE1))
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
  initial_conditions3 <- c(D0 = solution2[, "D02"][[2]],
                           Dm1 = solution2[, "Dm1"][[2]],
                           E1 = solution2[, "E1"][[2]])

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
  L1 <- solution3[, "D0"][[2]]
  logL1b <- log(L1)

  return(logL1b)

}





i <- 4
pars2 <- c(100, 11,0,2)
datalist[[i]]$branching_times <- c(4, 3, 1)

DAISIE_DE_logpES_max_min_age_coltime(datalist,
                                     i,
                                     pars1,
                                     methode = "lsodes",
                                     reltolint = 1e-16,
                                     abstolint = 1e-16)

DAISIE:::DAISIE_loglik_CS_choice(pars1 = pars1,
                                 pars2 = pars2,
                                 datalist = datalist,
                                 brts = datalist[[i]]$branching_times,
                                 stac = 9,
                                 missnumspec = datalist[[i]]$missing_species,
                                 methode = "odeint::runge_kutta_fehlberg78")

