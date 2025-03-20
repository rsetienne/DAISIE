#rm(list=ls())
library(deSolve)
library(DAISIE)
library(pracma)
###############################################################################
### fonction to calculate the likelihood of observing a non- endemic lineage
### with an unknown colonization time
###############################################################################

### Using D-E approach


# pars1[1] corresponds to the Cladogenesis rate
# pars1[2] corresponds to the Extinction rate of endemic lineages
# pars1[3] corresponds to the Extinction rate of non-endemic lineages
# pars1[4] = corresponds to the Colonization rate
# pars1[5] = corresponds to the Anagenesis rate
# if equal_extinction = TRUE, the extinction rates of endemic and non-endemic species are equal.
# else, the are estimated separately




DAISIE_DE_logpNE_unknown_coltime <- function(datalist,
                                          i,
                                          pars1,
                                          methode,
                                          rtol, 
                                          atol,
                                          equal_extinction = FALSE) {
  
  if (equal_extinction) {
    pars1[3] <- pars1[2]
  }
  
  
  t0 <- datalist[[1]]$island_age
  tp <- 0
  parameters <- pars1
  
  # Define system of equations for interval [t0, tp]
  interval1 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dD0 <- -pars1[4] * D0 + pars1[4] * Dm
      dDm <- -(pars1[5] + pars1[1] + pars1[3]) * Dm + (pars1[5] * E1 + pars1[1] * E1^2 + pars1[3]) * D0
      dE1 <- pars1[2] - (pars1[1] + pars1[2]) * E1 + pars1[1] * E1^2
      list(c(dD0, dDm, dE1))
    })
  }
  
  # Set initial conditions
  initial_conditions1 <- c(D0 = 0, Dm = 1, E1 = 0)
  
  # Time sequence for interval [t0, tp]
  time1 <- c(tp, t0)
  
  # Solve the system for interval [t0, tp]
  solution1 <- ode(y = initial_conditions1,
                   times = time1,
                   func = interval1,
                   parms = parameters,
                   method = methode,
                   rtol = rtol, 
                   atol = atol)
  
  # Extract log-likelihood
  L0 <- solution1[, "D0"][[2]]
  logL0b <- log(L0)
  return(logL0b)
}
