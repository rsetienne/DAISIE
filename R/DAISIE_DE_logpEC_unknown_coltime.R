#rm(list=ls())
library(deSolve)
library(DAISIE)
library(pracma)
###############################################################################
### fonction to calculate the likelihood of observing an endemic singleton lineage
### with the max age colonization time t1
###############################################################################

### Using D-E approach


# pars1[1] corresponds to the Cladogenesis rate
# pars1[2] corresponds to the Extinction rate of endemic lineages
# pars1[3] corresponds to the Extinction rate of non-endemic lineages
# pars1[4] = corresponds to the Colonization rate
# pars1[5] = corresponds to the Anagenesis rate
# if equal_extinction = TRUE, the extinction rates of endemic and non-endemic species are equal.
# else, the are estimated separately

### Using D-E approach
DAISIE_DE_logpEC_unknown_coltime <- function(datalist,
                                          i,
                                          pars1,
                                          methode,
                                          rtol, 
                                          atol,
                                          equal_extinction = FALSE) {
  if (equal_extinction) {
    pars1[3] <- pars1[2]
  }
  
  parameters <- pars1
  
  t0 <- datalist[[i]]$branching_times[1]
  t1 <- datalist[[i]]$branching_times[2]
  t2 <- datalist[[i]]$branching_times[3]
  tp <- 0
  ti <- sort(datalist[[i]]$branching_times)
  ti <- ti[1:(length(ti)-2)]
  parameters <- pars1
  
  # Define system of equations for interval [t2, tp]
  interval1 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dD1 <- -(pars1[1] + pars1[2]) * D1 + 2 * pars1[1] * D1 * E1
      dD0 <- -pars1[4] * D0 + pars1[4] * Dm
      dDm <- -(pars1[5] + pars1[1] + pars1[2]) * Dm + (pars1[5] * E1 + pars1[1] * E1^2 + pars1[2]) * D0 
      dE1 <- pars1[2] - (pars1[1] + pars1[2]) * E1 + pars1[1] * E1^2
      list(c(dD1, dD0, dDm, dE1))
    })
  }
  
  # Define system of equations for interval [t1, t2]
  interval2 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dD1 <- -(pars1[1] + pars1[2] ) * D1 + 2 * pars1[1] * D1 * E1
      dD0m <- -pars1[4] * D0m + pars1[4] * Dm
      dD0M <- -pars1[4] * D0M + pars1[4]*DM
      dDm <- -(pars1[5] + pars1[1] + pars1[2]) * Dm + (pars1[5] * E1 + pars1[1] * E1^2 + pars1[2]) * D0m 
      dDM <- -(pars1[5] + pars1[1] + pars1[2]) * DM + (pars1[5] * E1 + pars1[1] * E1^2 + pars1[2])* D0M + 
        (pars1[5] * D1 + 2 * pars1[1] * D1 * E1 ) * D0m 
      dE1 <- pars1[2] - (pars1[1] + pars1[2]) * E1 + pars1[1] * E1^2
      list(c(dD1, dD0m, dD0M, dDm, dDM, dE1))
    })
  }
  
  # Initial conditions
  number_of_species <- length(datalist[[i]]$branching_times) -1 
  number_of_missing_species <- datalist[[i]]$missing_species
  ro <- number_of_species / (number_of_missing_species + number_of_species) 

  initial_conditions1 <- c(D1 = ro, D0 = 1, Dm = 0, E1 = 1 - ro)
  
  solution0 <- ode(y = initial_conditions1,
                   times = c(0, ti),
                   func = interval1,
                   parms = parameters,
                   method = "lsodes",
                   rtol, atol)
  
  # Time sequences for interval [t2, tp]
  times <- rbind(c(0, ti[1:(length(ti) - 1)]), ti)
  
  for (idx in 1:length(ti)) {
    # Time sequence idx in interval [t2, tp]
    time1 <- times[, idx]
    
    # Solve the system for interval [t2, tp]
    solution1 <- ode(y = initial_conditions1,
                     times = time1,
                     func = interval1,
                     parms = parameters,
                     method = "lsodes",
                     rtol, atol)
    
    # Initial conditions
    initial_conditions1 <- c(D1 = pars1[1] * solution0[, "D1"][idx + 1] * solution1[, "D1"][2],
                             D0 = 1, Dm = 0, E1 = solution0[, "E1"][idx + 1])
  }
  
  # Initial conditions
  initial_conditions2 <- c(D1 = initial_conditions1["D1"][[1]],
                           D0m = solution0[, "D0"][length(ti) + 1],
                           D0M = 0,
                           Dm = solution0[, "Dm"][length(ti) + 1],
                           DM = initial_conditions1["D1"][[1]] * solution0[, "D0"][length(ti)+1],
                           E1 = initial_conditions1["E1"][[1]])
  
  # Time sequence for interval [t1, t2]
  time2 <- c(t2, t1)
  
  # Solve the system for interval [t2, tp]
  solution2 <- ode(y = initial_conditions2,
                   times = time2,
                   func = interval2,
                   parms = parameters,
                   method = "lsodes",
                   rtol, atol)
  
  
  # Extract log-likelihood
  Lk <- (solution2[, "D0M"][[2]])
  logLkb <- log(Lk)
  return(logLkb)
}
i <- 23
datalist[[i]]$branching_times <- c(52, 51.99999, 45)

DAISIE_DE_logpEC_unknown_coltime(datalist,i,pars1,
                                          methode = "lsodes",
                                          rtol = 1e-16, 
                                          atol = 1e-16,
                                          equal_extinction = FALSE)

DAISIE_DE_logpEC_max_age_coltime(datalist,i,pars1,
                                             methode = "lsodes",
                                             rtol = 1e-16, 
                                             atol = 1e-16,
                                             equal_extinction = TRUE)


DAISIE:::DAISIE_loglik_all(pars11, pars2, datalist, 
                             methode = "odeint::runge_kutta_fehlberg78",
                             abstolint = 1e-16, reltolint = 1e-16)
