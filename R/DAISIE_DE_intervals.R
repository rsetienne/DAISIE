# pars1:
# 1: lambda_c
# 2: mu (endemic)
# 3: mu2 (non-endemic)
# 4: gamma
# 5: lambda_a


# interval functions
# DAISIE_DE_logpNE_max_age_coltime
# DAISIE_DE_logpNE_max_min_age_coltime
#' @keywords internal
interval3_NE <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    lambdac <- parameters[1]
    mu      <- parameters[2]
    gamma   <- parameters[4]
    lambdaa <- parameters[5]

    dDM1 <- -(lambdac + mu + lambdaa + gamma) * DM1 + (mu + lambdaa * E + lambdac * E * E) * DA2 + gamma * DM2
    dDM2 <- -(lambdac + mu + lambdaa) * DM2 + (mu + lambdaa * E + lambdac * E * E) * DA2

    dE <- mu - (mu + lambdac) * E + lambdac * E * E
    dDA2 <- -gamma * DA2 + gamma * DM2

    return(list(c(dDM1, dDM2, dE, dDA2)))
  })
}

# Define system of equations for interval [t2, tp]
# DAISIE_DE_logpEC_mainland
# DAISIE_DE_logpEC_max_age_coltime_and_mainland
# DAISIE_DE_logpEC_max_age_coltime
# DAISIE_DE_logpEC
#' @keywords internal
interval2_EC <- function(t, state, pars1) {
  with(as.list(c(state, pars1)), {
    lambdac <- parameters[1]
    mu      <- parameters[2]
    gamma   <- parameters[4]
    lambdaa <- parameters[5]

    dDE  <- -(lambdac + mu) * DE + 2 * lambdac * DE * E
    dDM3 <- -(lambdac + mu + lambdaa) * DM3 + (mu + lambdaa * E + lambdac * E * E) * DA3
    dE   <- mu - (mu + lambdac) * E + lambdac * E * E
    dDA3 <- -gamma * DA3 + gamma * DM3

    return(list(c(dDE, dDM3, dE, dDA3)))
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
interval4 <- function(t, state, pars1) {
  with(as.list(c(state, pars1)), {
    lambdac <- parameters[1]
    mu      <- parameters[2]
    gamma   <- parameters[4]
    lambdaa <- parameters[5]

    dDA1 <- -gamma * DA1 + gamma * DM1

    dDM1 <- -(lambdac + mu + lambdaa) * DM1 + (mu + lambdaa * E + lambdac * E * E) * DA1
    dE <- mu - (mu + lambdac) * E + lambdac * E * E
    return(list(c( dDA1, dDM1, dE)))
  })
}


# Define system of equations for interval [t1, t2]
# DAISIE_DE_logpEC_max_age_coltime_and_mainland
# DAISIE_DE_logpEC_max_age_coltime
# DAISIE_DE_logpES_max_age_coltime_and_mainland
# DAISIE_DE_logpES_max_age_coltime
# DAISIE_DE_logpES_max_min_age_coltime
#' @keywords internal
interval3_ES <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    lambdac <- parameters[1]
    mu      <- parameters[2]
    gamma   <- parameters[4]
    lambdaa <- parameters[5]


    dDE <- -(lambdac + mu) * DE + 2 * lambdac * DE * E

    dDM1 <- -(lambdac + mu + lambdaa + gamma) * DM1 + gamma * DM2 + (mu + lambdaa * E + lambdac * E * E) * DA2
    dDM2 <- -(lambdac + mu + lambdaa) * DM2 + (mu + lambdaa * E + lambdac * E * E) * DA2 + (lambdaa * DE + 2 * lambdac * DE * E) * DA3
    dDM3 <- -(lambdac + mu + lambdaa) * DM3 + (mu + lambdaa * E + lambdac * E * E) * DA3

    dE <- mu - (mu + lambdac) * E + lambdac * E * E

    dDA2 <- -gamma * DA2 + gamma * DM2
    dDA3 <- -gamma * DA3 + gamma * DM3

    return(list(c(dDE, dDM1, dDM2, dDM3, dE, dDA2, dDA3)))
  })
}

# Define system of equations for interval [t1, tp]
# DAISIE_DE_logpES_mainland
#' @keywords internal
interval2_ES <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    lambdac <- parameters[1]
    mu      <- parameters[2]
    gamma   <- parameters[4]
    lambdaa <- parameters[5]

    dDE <- -(lambdac + mu) * DE + 2 * lambdac * DE * E

    dDM2 <- -(lambdac + mu + gamma + lambdaa) * DM2 + (lambdaa * DE + 2 * lambdac * DE * E) * DA3

    dDM3 <- -(lambdac + mu + lambdaa) * DM3 + (mu + lambdaa * E + lambdac * E * E) * DA3

    dE <- mu - (mu + lambdac) * E + lambdac * E * E

    dDA3 <- -gamma * DA3 + gamma * DM3


    return(list(c(dDE, dDM2, dDM3, dE, dDA3)))

  })
}


# Define system of equations for interval [t1, tp]
# DAISIE_DE_logpNE_max_min_age_coltime
# DAISIE_DE_logpNE
#' @keywords internal
interval2_NE <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    lambdac <- parameters[1]
    mu      <- parameters[2]
    gamma   <- parameters[4]
    lambdaa <- parameters[5]


    dDM2 <- -(lambdac + mu + gamma + lambdaa) * DM2

    dE <- mu - (mu + lambdac) * E + lambdac * E * E

    return(list(c(dDM2, dE)))
  })
}
