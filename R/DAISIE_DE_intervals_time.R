# pars1:
# 1: lambda_c
# 2: mu_E (endemic)
# 3: mu_NE (non-endemic)
# 4: gamma
# 5: lambda_a


# interval functions
# DAISIE_DE_logpNE_max_age_coltime
# DAISIE_DE_logpNE_max_min_age_coltime
#' @keywords internal
#'

island_area_vector <- function(timeval,
                               area_pars,
                               island_ontogeny,
                               sea_level,
                               total_time,
                               peak) {
  # Constant
  if (island_ontogeny == 0 || is.na(island_ontogeny)) {
    if (area_pars[1] != 1 || is.null(area_pars[1])) {
      warning("Constant ontogeny requires a maximum area of 1.")
    }
    return(1)
  } else { # Ontogeny
    area_pars <- create_area_pars(max_area = area_pars[1],
                                  current_area = area_pars[2],
                                  proportional_peak_t = area_pars[3],
                                  total_island_age = area_pars[4],
                                  sea_level_amplitude = area_pars[5],
                                  sea_level_frequency = area_pars[6],
                                  island_gradient_angle = area_pars[7])



    area <- island_area(
      timeval = timeval,
      area_pars = area_pars,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      total_time = total_time,
      peak = peak
    )
    return(area)
  }
}

interval3_NE_time_dep <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    lambdac_0 <- parameters[1]
    mu_E_0    <- parameters[2]
    mu_NE_0   <- parameters[3]
    gamma_0   <- parameters[4]
    lambdaa_0 <- parameters[5]
    d <- parameters[6]
    x_hyperpar <- parameters[7]
    area_pars <- parameters[8:14]
    island_ontogeny <- parameters[15]
    sea_level <- parameters[16]
    total_time <- parameters[17]
    peak <- parameters[18]
    kk <- parameters[19]
    ddep <- parameters[20]



    area <- island_area_vector(

      timeval = total_time - t,
      area_pars = area_pars,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      total_time = total_time,
      peak = peak
    )

    lambdac <- get_clado_rate_per_capita(
      lac = lambdac_0,
      d = d,
      A = area,
      K = Inf,
      num_spec = 0
    )
    mu_E <- get_ext_rate_per_capita(
      mu = mu_E_0,
      x = x_hyperpar,
      A = area,
      extcutoff = 1000000
    )

    mu_NE <- get_ext_rate_per_capita(
      mu = mu_NE_0,
      x = x_hyperpar,
      A = area,
      extcutoff = 1000000
    )

    gamma <- get_immig_rate_per_capita(
      gam = gamma_0,
      A = area,
      num_spec = 0,
      K = Inf
    )
    lambdaa <- lambdaa_0


    dDM1 <- -(lambdac + mu_NE + lambdaa + gamma) * DM1 + (mu_NE + lambdaa * E + lambdac * E * E) * DA2 + gamma * DM2
    dDM2 <- -(lambdac + mu_NE + lambdaa) * DM2 + (mu_NE + lambdaa * E + lambdac * E * E) * DA2

    dE <- mu_E - (mu_E + lambdac) * E + lambdac * E * E
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
interval2_EC_time_dep <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    lambdac_0 <- parameters[1]
    mu_E_0    <- parameters[2]
    mu_NE_0   <- parameters[3]
    gamma_0   <- parameters[4]
    lambdaa_0 <- parameters[5]
    d <- parameters[6]
    x_hyperpar <- parameters[7]
    area_pars <- parameters[8:14]
    island_ontogeny <- parameters[15]
    sea_level <- parameters[16]
    total_time <- parameters[17]
    peak <- parameters[18]
    kk <- parameters[19]
    ddep <- parameters[20]



    area <- island_area_vector(

      timeval = total_time - t,
      area_pars = area_pars,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      total_time = total_time,
      peak = peak
    )

    lambdac <- get_clado_rate_per_capita(
      lac = lambdac_0,
      d = d,
      A = area,
      K = Inf,
      num_spec = 0
    )
    mu_E <- get_ext_rate_per_capita(
      mu = mu_E_0,
      x = x_hyperpar,
      A = area,
      extcutoff = 1000000
    )

    mu_NE <- get_ext_rate_per_capita(
      mu = mu_NE_0,
      x = x_hyperpar,
      A = area,
      extcutoff = 1000000
    )

    gamma <- get_immig_rate_per_capita(
      gam = gamma_0,
      A = area,
      num_spec = 0,
      K = Inf
    )
    lambdaa <- lambdaa_0


    dDE  <- -(lambdac + mu_E) * DE + 2 * lambdac * DE * E
    dDM3 <- -(lambdac + mu_NE + lambdaa) * DM3 + (mu_NE + lambdaa * E + lambdac * E * E) * DA3
    dE   <- mu_E - (mu_E + lambdac) * E + lambdac * E * E
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
interval4_time_dep <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    lambdac_0 <- parameters[1]
    mu_E_0    <- parameters[2]
    mu_NE_0   <- parameters[3]
    gamma_0   <- parameters[4]
    lambdaa_0 <- parameters[5]
    d <- parameters[6]
    x_hyperpar <- parameters[7]
    area_pars <- parameters[8:14]
    island_ontogeny <- parameters[15]
    sea_level <- parameters[16]
    total_time <- parameters[17]
    peak <- parameters[18]
    kk <- parameters[19]
    ddep <- parameters[20]



    area <- island_area_vector(

      timeval = total_time - t,
      area_pars = area_pars,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      total_time = total_time,
      peak = peak
    )

    lambdac <- get_clado_rate_per_capita(
      lac = lambdac_0,
      d = d,
      A = area,
      K = Inf,
      num_spec = 0
    )
    mu_E <- get_ext_rate_per_capita(
      mu = mu_E_0,
      x = x_hyperpar,
      A = area,
      extcutoff = 1000000
    )

    mu_NE <- get_ext_rate_per_capita(
      mu = mu_NE_0,
      x = x_hyperpar,
      A = area,
      extcutoff = 1000000
    )

    gamma <- get_immig_rate_per_capita(
      gam = gamma_0,
      A = area,
      num_spec = 0,
      K = Inf
    )
    lambdaa <- lambdaa_0


    dDA1 <- -gamma * DA1 + gamma * DM1

    dDM1 <- -(lambdac + mu_NE + lambdaa) * DM1 + (mu_NE + lambdaa * E + lambdac * E * E) * DA1
    dE <- mu_E - (mu_E + lambdac) * E + lambdac * E * E
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
interval3_ES_time_dep <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    lambdac_0 <- parameters[1]
    mu_E_0    <- parameters[2]
    mu_NE_0   <- parameters[3]
    gamma_0   <- parameters[4]
    lambdaa_0 <- parameters[5]
    d <- parameters[6]
    x_hyperpar <- parameters[7]
    area_pars <- parameters[8:14]
    island_ontogeny <- parameters[15]
    sea_level <- parameters[16]
    total_time <- parameters[17]
    peak <- parameters[18]
    kk <- parameters[19]
    ddep <- parameters[20]



    area <- island_area_vector(

      timeval = total_time - t,
      area_pars = area_pars,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      total_time = total_time,
      peak = peak
    )

    lambdac <- get_clado_rate_per_capita(
      lac = lambdac_0,
      d = d,
      A = area,
      K = Inf,
      num_spec = 0
    )
    mu_E <- get_ext_rate_per_capita(
      mu = mu_E_0,
      x = x_hyperpar,
      A = area,
      extcutoff = 1000000
    )

    mu_NE <- get_ext_rate_per_capita(
      mu = mu_NE_0,
      x = x_hyperpar,
      A = area,
      extcutoff = 1000000
    )

    gamma <- get_immig_rate_per_capita(
      gam = gamma_0,
      A = area,
      num_spec = 0,
      K = Inf
    )
    lambdaa <- lambdaa_0



    dDE <- -(lambdac + mu_E) * DE + 2 * lambdac * DE * E

    dDM1 <- -(lambdac + mu_NE + lambdaa + gamma) * DM1 + gamma * DM2 + (mu_NE + lambdaa * E + lambdac * E * E) * DA2
    dDM2 <- -(lambdac + mu_NE + lambdaa) * DM2 + (mu_NE + lambdaa * E + lambdac * E * E) * DA2 + (lambdaa * DE + 2 * lambdac * DE * E) * DA3
    dDM3 <- -(lambdac + mu_NE + lambdaa) * DM3 + (mu_NE + lambdaa * E + lambdac * E * E) * DA3

    dE <- mu_E - (mu_E + lambdac) * E + lambdac * E * E

    dDA2 <- -gamma * DA2 + gamma * DM2
    dDA3 <- -gamma * DA3 + gamma * DM3

    return(list(c(dDE, dDM1, dDM2, dDM3, dE, dDA2, dDA3)))
  })
}

# Define system of equations for interval [t1, tp]
# DAISIE_DE_logpES_mainland
#' @keywords internal
interval2_ES_time_dep <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    lambdac_0 <- parameters[1]
    mu_E_0    <- parameters[2]
    mu_NE_0   <- parameters[3]
    gamma_0   <- parameters[4]
    lambdaa_0 <- parameters[5]
    d <- parameters[6]
    x_hyperpar <- parameters[7]
    area_pars <- parameters[8:14]
    island_ontogeny <- parameters[15]
    sea_level <- parameters[16]
    total_time <- parameters[17]
    peak <- parameters[18]
    kk <- parameters[19]
    ddep <- parameters[20]



    area <- island_area_vector(

      timeval = total_time - t,
      area_pars = area_pars,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      total_time = total_time,
      peak = peak
    )

    lambdac <- get_clado_rate_per_capita(
      lac = lambdac_0,
      d = d,
      A = area,
      K = Inf,
      num_spec = 0
    )
    mu_E <- get_ext_rate_per_capita(
      mu = mu_E_0,
      x = x_hyperpar,
      A = area,
      extcutoff = 1000000
    )

    mu_NE <- get_ext_rate_per_capita(
      mu = mu_NE_0,
      x = x_hyperpar,
      A = area,
      extcutoff = 1000000
    )

    gamma <- get_immig_rate_per_capita(
      gam = gamma_0,
      A = area,
      num_spec = 0,
      K = Inf
    )
    lambdaa <- lambdaa_0


    dDE <- -(lambdac + mu_E) * DE + 2 * lambdac * DE * E

    dDM2 <- -(lambdac + mu_NE + gamma + lambdaa) * DM2 + (lambdaa * DE + 2 * lambdac * DE * E) * DA3

    dDM3 <- -(lambdac + mu_NE + lambdaa) * DM3 + (mu_NE + lambdaa * E + lambdac * E * E) * DA3

    dE <- mu_E - (mu_E + lambdac) * E + lambdac * E * E

    dDA3 <- -gamma * DA3 + gamma * DM3


    return(list(c(dDE, dDM2, dDM3, dE, dDA3)))

  })
}


# Define system of equations for interval [t1, tp]
# DAISIE_DE_logpNE_max_min_age_coltime
# DAISIE_DE_logpNE
#' @keywords internal
interval2_NE_time_dep <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    lambdac_0 <- parameters[1]
    mu_E_0    <- parameters[2]
    mu_NE_0   <- parameters[3]
    gamma_0   <- parameters[4]
    lambdaa_0 <- parameters[5]
    d <- parameters[6]
    x_hyperpar <- parameters[7]
    area_pars <- parameters[8:14]
    island_ontogeny <- parameters[15]
    sea_level <- parameters[16]
    total_time <- parameters[17]
    peak <- parameters[18]
    kk <- parameters[19]
    ddep <- parameters[20]



    area <- island_area_vector(

      timeval = total_time - t,
      area_pars = area_pars,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      total_time = total_time,
      peak = peak
    )

    lambdac <- get_clado_rate_per_capita(
      lac = lambdac_0,
      d = d,
      A = area,
      K = Inf,
      num_spec = 0
    )
    mu_E <- get_ext_rate_per_capita(
      mu = mu_E_0,
      x = x_hyperpar,
      A = area,
      extcutoff = 1000000
    )

    mu_NE <- get_ext_rate_per_capita(
      mu = mu_NE_0,
      x = x_hyperpar,
      A = area,
      extcutoff = 1000000
    )

    gamma <- get_immig_rate_per_capita(
      gam = gamma_0,
      A = area,
      num_spec = 0,
      K = Inf
    )
    lambdaa <- lambdaa_0



    dDM2 <- -(lambdac + mu_NE + gamma + lambdaa) * DM2

    dE <- mu_E - (mu_E + lambdac) * E + lambdac * E * E

    return(list(c(dDM2, dE)))
  })
}
