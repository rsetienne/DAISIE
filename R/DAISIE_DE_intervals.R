# interval functions
# DAISIE_DE_logpNE_max_age_coltime
# DAISIE_DE_logpNE_max_min_age_coltime
#' @keywords internal
interval1 <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dDA <-  -pars1[4] * DA + pars1[4] * Dm2
    dDm1 <- -(pars1[5] + pars1[1] + pars1[3] + pars1[4]) * Dm1 +
      (pars1[5] * E + pars1[1] * E^2 + pars1[3]) * DA + pars1[4] * Dm2
    dDm2 <- -(pars1[5] + pars1[1] + pars1[3]) * Dm2 +
      (pars1[5] * E + pars1[1] * E^2 + pars1[3]) * DA
    dE <-  pars1[2] - (pars1[1] + pars1[2]) * E + pars1[1] * E^2
    list(c(dDA, dDm1, dDm2, dE))
  })
}

# Define system of equations for interval [t2, tp]
# DAISIE_DE_logpEC_mainland
# DAISIE_DE_logpEC_max_age_coltime_and_mainland
# DAISIE_DE_logpEC_max_age_coltime
# DAISIE_DE_logpEC
#' @keywords internal
interval1_3 <- function(t, state, pars1) {
  with(as.list(c(state, pars1)), {
    dDE <- -(pars1[1] + pars1[2]) * DE + 2 * pars1[1] * DE * E
    dDA3 <- -pars1[4] * DA3 + pars1[4] * Dm3
    dDm3 <- -(pars1[5] + pars1[1] + pars1[3]) * Dm3 + (pars1[5] * E + pars1[1] * E^2 + pars1[3]) * DA3
    dE <- pars1[2] - (pars1[1] + pars1[2]) * E + pars1[1] * E^2
    list(c(dDE, dDA3, dDm3, dE))
  })
}

# Define system of equations for interval [t1, tp]
# DAISIE_DE_logpES
# DAISIE_DE_logpEC
# DAISIE_DE_logpEC_mainland <<-- checkecheckehcek!
# DAISIE_DE_logpES_max_min_age_coltime
#' @keywords internal
interval2 <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dDE <- -(pars1[1] + pars1[2]) * DE + 2 * pars1[1] * DE * E
    dDA3 <- -pars1[4] * DA3 + pars1[4] * Dm3
    dDm3 <- -(pars1[5] + pars1[1] + pars1[3]) * Dm3 + (pars1[5] * E + pars1[1] * E^2 + pars1[3]) * DA3
    dDm2 <- -(pars1[5] + pars1[1] + pars1[3] + pars1[4]) * Dm2 + (pars1[5] * DE + 2 * pars1[1] * DE * E) * DA3
    dE <- pars1[2] - (pars1[1] + pars1[2]) * E + pars1[1] * E^2
    list(c(dDE, dDA3, dDm3, dDm2, dE))
  })
}

# Define system of equations for interval [t0, t1]
# DAISIE_DE_logpEC_mainland
# DAISIE_DE_logpEC_max_age_coltime_and_mainland
# DAISIE_DE_logpEC_max_age_coltime
# DAISIE_DE_logpEC
# DAISIE_DE_logpES_mainland
# DAISIE_DE_logpES_max_age_coltime_and_mainland
# DAISIE_DE_logpES_max_age_coltime
# DAISIE_DE_logpES_max_min_age_coltime
# DAISIE_DE_logpNE_unknown_coltime
# DAISIE_DE_logpNE
# DAISIE_DE_logp0
# DAISIE_DE_logpES
# DAISIE_DE_logpNE_max_age_coltime
# DAISIE_DE_logpNE_max_min_age_coltime
#' @keywords internal
interval3 <- function(t, state, pars1) {
  with(as.list(c(state, pars1)), {
    dDA1 <- -pars1[4] * DA1 + pars1[4] * Dm1
    dDm1 <- -(pars1[5] + pars1[1] + pars1[3]) * Dm1 + (pars1[5] * E + pars1[1] * E^2 + pars1[3]) * DA1
    dE <- pars1[2] - (pars1[1] + pars1[2]) * E + pars1[1] * E^2
    list(c(dDA1, dDm1, dE))
  })
}


# Define system of equations for interval [t1, t2]
# DAISIE_DE_logpEC_max_age_coltime_and_mainland
# DAISIE_DE_logpEC_max_age_coltime
# DAISIE_DE_logpES_max_age_coltime_and_mainland
# DAISIE_DE_logpES_max_age_coltime
# DAISIE_DE_logpES_max_min_age_coltime
#' @keywords internal
interval2_4 <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dDE <- -(pars1[1] + pars1[2]) * DE + 2 * pars1[1] * DE * E
    dDA2 <- -pars1[4] * DA2 + pars1[4] * Dm2
    dDA3 <- -pars1[4] * DA3 + pars1[4] * Dm3
    dDm1 <- -(pars1[5] + pars1[1] + pars1[3] + pars1[4]) * Dm1 +
      (pars1[3] + pars1[5] * E + pars1[1] * E^2) * DA2 + pars1[4] * Dm2
    dDm2 <- -(pars1[5] + pars1[1] + pars1[3]) * Dm2 + (pars1[3] + pars1[5] * E + pars1[1] * E^2) * DA2 +
      (pars1[5] * DE + 2 * pars1[1] * DE * E ) * DA3
    dDm3 <- -(pars1[5] + pars1[1] + pars1[3]) * Dm3 + (pars1[3] + pars1[5] * E + pars1[1] * E^2) * DA3
    dE <- pars1[2] - (pars1[1] + pars1[2]) * E + pars1[1] * E^2
    list(c(dDE, dDA2, dDA3, dDm1, dDm2, dDm3, dE))
  })
}

# Define system of equations for interval [t2, tp]
# DAISIE_DE_logpEC_unknown_coltime
#' @keywords internal
interval1_6 <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dD1 <- -(pars1[1] + pars1[2]) * D1 + 2 * pars1[1] * D1 * E1
    dD0 <- -pars1[4] * D0 + pars1[4] * Dm
    dDm <- -(pars1[5] + pars1[1] + pars1[2]) * Dm + (pars1[5] * E1 + pars1[1] * E1^2 + pars1[2]) * D0
    dE1 <- pars1[2] - (pars1[1] + pars1[2]) * E1 + pars1[1] * E1^2
    list(c(dD1, dD0, dDm, dE1))
  })
}

# Define system of equations for interval [t1, t2]
# DAISIE_DE_logpEC_unknown_coltime
#' @keywords internal
interval2_6 <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dD1 <- -(pars1[1] + pars1[2] ) * D1 + 2 * pars1[1] * D1 * E1
    dD0m <- -pars1[4] * D0m + pars1[4] * Dm
    dD0M <- -pars1[4] * D0M + pars1[4] * DM
    dDm <- -(pars1[5] + pars1[1] + pars1[2]) * Dm + (pars1[5] * E1 + pars1[1] * E1^2 + pars1[2]) * D0m
    dDM <- -(pars1[5] + pars1[1] + pars1[2]) * DM + (pars1[5] * E1 + pars1[1] * E1^2 + pars1[2]) * D0M + (pars1[5] * D1 + 2 * pars1[1] * D1 * E1 ) * D0m
    dE1 <- pars1[2] - (pars1[1] + pars1[2]) * E1 + pars1[1] * E1^2
    list(c(dD1, dD0m, dD0M, dDm, dDM, dE1))
  })
}


# Define system of equations for interval [t1, tp]
# DAISIE_DE_logpES_unknown_coltime
#' @keywords internal
interval1_12 <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dDE <- -(pars1[1] + pars1[2]) * DE + 2 * pars1[1] * DE * E
    dDA3 <- -pars1[4] * DA3 + pars1[4] * Dm3
    dDA2 <- -pars1[4] * DA2 + pars1[4] * Dm2
    dDm3 <- -(pars1[5] + pars1[1] + pars1[3]) * Dm3 + (pars1[5] * E + pars1[1] * E^2 + pars1[3]) * DA3
    dDm2 <- -(pars1[5] + pars1[1] + pars1[3]) * Dm2 + (pars1[5] * DE + 2 * pars1[1] * DE * E) * DA3 + (pars1[3] + pars1[5] * E + pars1[1] * E^2) * DA2
    dE <- pars1[2] - (pars1[1] + pars1[2]) * E + pars1[1] * E^2
    list(c(dDE, dDA3, dDA2, dDm3, dDm2, dE))
  })
}

# Define system of equations for interval [t1, tp]
# DAISIE_DE_logpES_mainland
#' @keywords internal
interval1_8 <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dDE <- -(pars1[1] + pars1[2]) * DE + 2 * pars1[1] * DE * E
    dDA3 <- -pars1[4] * DA3 + pars1[4] * Dm
    dDm3 <- -(pars1[5] + pars1[1] + pars1[3]) * Dm + (pars1[5] * E + pars1[1] * E^2 + pars1[3]) * DA3
    dDm2 <- -(pars1[5] + pars1[1] + pars1[3] + pars1[4]) * Dm2 + (pars1[5] * D1 + 2 * pars1[1] * DE * E) * DA3
    dE <- pars1[2] - (pars1[1] + pars1[2]) * E + pars1[1] * E^2
    list(c(dDE, dDA3, dDm3, dDm2, dE))
  })
}


# Define system of equations for interval [t1, tp]
# DAISIE_DE_logpNE_max_min_age_coltime
# DAISIE_DE_logpNE
#' @keywords internal
interval1_13 <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dDm2 <- -(pars1[5] + pars1[1] + pars1[3] + pars1[4]) * Dm2
    dE <-  pars1[2] - (pars1[1] + pars1[2]) * E + pars1[1] * E^2
    list(c(dDm2, dE))
  })
}
