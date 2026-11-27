# Diversity-dependent versions of DAISIE_DE_logp0, DAISIE_DE_logpNE,
# DAISIE_DE_logpES and DAISIE_DE_logpEC. They follow the same intervals as the
# DAISIE_DE functions; the equations are built in DAISIE_DE_DD_intervals.R.
#
# pars1 = c(lambda_c, mu, K, gamma, lambda_a), as in DAISIE_loglik_CS.
# pars2[1] = largest number of hidden species kept when the clade size is not
#            limited by K (as lx in DAISIE; with ddep = 11 the limit
#            ceiling(K) is used when it is smaller)
# pars2[2] = ddep (0, 1, 11, 2 or 21), as in DAISIE
# sampling = 'n': missnumspec is the exact number of missing species
#            'rho': each species is sampled with probability
#            S / (S + missnumspec), as in DAISIE_DE with sampling = 'rho'

# Hidden species at the present: a vector over n = 0, 1, ...
DAISIE_DE_DD_hidden_init <- function(missnumspec, S, sampling, n_cap) {
  if (sampling == 'rho' && S > 0 && missnumspec > 0) {
    rho <- S / (S + missnumspec)
    return((1 - rho)^(0:(n_cap + missnumspec)))
  }
  h <- rep(0, missnumspec + 1)
  h[missnumspec + 1] <- 1
  h
}

# Correction for the sampling of the S observed species
DAISIE_DE_DD_sampling_term <- function(missnumspec, S, sampling) {
  if (S == 0 || missnumspec == 0) return(0)
  if (sampling == 'rho') return(S * (log(S) - log(S + missnumspec)))
  # number of ways to choose which of the S + missnumspec species are missing
  return(-lchoose(S + missnumspec, missnumspec))
}

# From the colonization time t1 back to the island age t0 (interval4).
# DA1, DM1: vectors over n at t1. Returns log(DA1[n = 0]) at t0.
DAISIE_DE_DD_interval4_to_island_age <- function(DA1, DM1, t1, t0, rates,
                                                 n_cap, methode, atol, rtol) {
  system4 <- DAISIE_DE_DD_system("interval4", 0, rates, n_cap)
  y <- DAISIE_DE_DD_state(system4, list(DA1 = DA1, DM1 = DM1))
  y <- DAISIE_DE_DD_solve_branch(system4, y, c(t1, t0), methode, atol, rtol)
  log(DAISIE_DE_DD_get(system4, y, "DA1")[1])
}

# Colonization at the (known) colonization time: DM2 is the colonist.
# The mainland species was absent (new colonization) or present as a
# non-endemic species (recolonization).
DAISIE_DE_DD_colonization <- function(DM2, rates) {
  n <- seq_along(DM2) - 1
  list(DA1 = rates$colon(n) * DM2,
       DM1 = rates$recol(n + 1) * DM2)
}

DAISIE_DE_DD_logp0 <- function(island_age,
                               pars1,
                               pars2 = c(100, 11, 0, 0),
                               methode = "odeint::runge_kutta_cash_karp54",
                               reltolint = 1e-10,
                               abstolint = 1e-16,
                               dd_recolonization = FALSE,
                               mu_NE = pars1[2]) {
  rates <- DAISIE_DE_DD_rates(pars1, pars2[2], dd_recolonization, mu_NE)
  # nothing on the island at present: no hidden species, no non-endemic
  DAISIE_DE_DD_interval4_to_island_age(DA1 = 1, DM1 = 0, t1 = 0, t0 = island_age,
                                       rates = rates, n_cap = pars2[1],
                                       methode = methode,
                                       atol = abstolint, rtol = reltolint)
}

DAISIE_DE_DD_logpNE <- function(brts,
                                missnumspec = 0,
                                pars1,
                                pars2 = c(100, 11, 0, 0),
                                stac,
                                methode = "odeint::runge_kutta_cash_karp54",
                                reltolint = 1e-10,
                                abstolint = 1e-16,
                                dd_recolonization = FALSE,
                                mu_NE = pars1[2]) {
  if (!(stac %in% c(1, 4, 8))) {
    stop("NE only supports stac values of 1, 4 and 8")
  }
  t0 <- brts[1]
  t1 <- brts[2]
  t2 <- brts[3]
  n_cap <- pars2[1] + missnumspec
  rates <- DAISIE_DE_DD_rates(pars1, pars2[2], dd_recolonization, mu_NE)
  hidden <- DAISIE_DE_DD_hidden_init(missnumspec, 0, 'n', n_cap)

  if (stac == 4 || stac == 8) {
    # known colonization time (stac 4) or minimum colonization time t2 (stac 8):
    # the non-endemic species stays as it is from t1 (or t2) to the present
    t_end <- if (stac == 4) t1 else t2
    system2 <- DAISIE_DE_DD_system("interval2_NE", 0, rates, n_cap)
    y <- DAISIE_DE_DD_state(system2, list(DM2 = hidden))
    y <- DAISIE_DE_DD_solve_branch(system2, y, c(0, t_end), methode, abstolint, reltolint)
    DM2 <- DAISIE_DE_DD_get(system2, y, "DM2")
    if (stac == 4) {
      col <- DAISIE_DE_DD_colonization(DM2, rates)
      return(DAISIE_DE_DD_interval4_to_island_age(col$DA1, col$DM1, t1, t0, rates,
                                                  n_cap, methode, abstolint, reltolint))
    }
    initial <- list(DM2 = DM2)
    time3 <- c(t2, t1)
  } else {
    # maximum colonization time (stac 1)
    initial <- list(DM2 = hidden)
    time3 <- c(0, t1)
  }
  system3 <- DAISIE_DE_DD_system("interval3_NE", 0, rates, n_cap)
  y <- DAISIE_DE_DD_state(system3, initial)
  y <- DAISIE_DE_DD_solve_branch(system3, y, time3, methode, abstolint, reltolint)
  DAISIE_DE_DD_interval4_to_island_age(DAISIE_DE_DD_get(system3, y, "DA2"),
                                       DAISIE_DE_DD_get(system3, y, "DM1"),
                                       t1, t0, rates, n_cap, methode,
                                       abstolint, reltolint)
}

# From the stem of the observed clade (one lineage, k = 1) at time t_start back
# to the island age. DA3, DM3, DM2: vectors over n at t_start.
# stac 2, 3: known colonization time t1; stac 5, 6, 7, 9: maximum
# colonization time t1.
DAISIE_DE_DD_stem_to_island_age <- function(DA3, DM3, DM2, t_start, t1, t0,
                                            stac, rates, n_cap, methode,
                                            atol, rtol) {
  if (stac %in% c(2, 3)) {
    system2 <- DAISIE_DE_DD_system("interval2_ES", 1, rates, n_cap)
    y <- DAISIE_DE_DD_state(system2, list(DA3 = DA3, DM3 = DM3, DM2 = DM2))
    y <- DAISIE_DE_DD_solve_branch(system2, y, c(t_start, t1), methode, atol, rtol)
    col <- DAISIE_DE_DD_colonization(DAISIE_DE_DD_get(system2, y, "DM2"), rates)
    DA1 <- col$DA1
    DM1 <- col$DM1
  } else {
    system3 <- DAISIE_DE_DD_system("interval3_ES", 1, rates, n_cap)
    y <- DAISIE_DE_DD_state(system3, list(DA3 = DA3, DM3 = DM3, DM2 = DM2))
    y <- DAISIE_DE_DD_solve_branch(system3, y, c(t_start, t1), methode, atol, rtol)
    DA1 <- DAISIE_DE_DD_get(system3, y, "DA2")
    DM1 <- DAISIE_DE_DD_get(system3, y, "DM1")
  }
  DAISIE_DE_DD_interval4_to_island_age(DA1, DM1, t1, t0, rates, n_cap,
                                       methode, atol, rtol)
}

DAISIE_DE_DD_logpES <- function(brts,
                                missnumspec = 0,
                                pars1,
                                pars2 = c(100, 11, 0, 0),
                                stac,
                                methode = "odeint::runge_kutta_cash_karp54",
                                reltolint = 1e-10,
                                abstolint = 1e-16,
                                sampling = 'n',
                                dd_recolonization = FALSE,
                                mu_NE = pars1[2]) {
  if (!(stac %in% c(2, 3, 5, 7, 9))) {
    stop("stac must be 2, 3, 5, 7 or 9 for this function.")
  }
  t0 <- brts[1]
  t1 <- brts[2]
  t2 <- brts[3]
  S <- 1
  n_cap <- pars2[1] + missnumspec
  rates <- DAISIE_DE_DD_rates(pars1, pars2[2], dd_recolonization, mu_NE)
  hidden <- DAISIE_DE_DD_hidden_init(missnumspec, S, sampling, n_cap)
  # at present: the observed endemic species, with (stac 3, 7) or without
  # the non-endemic mainland species on the island
  if (stac %in% c(3, 7)) {
    initial <- list(DM3 = hidden)
  } else {
    initial <- list(DA3 = hidden)
  }
  t_start <- 0
  stac_stem <- stac
  if (stac == 9) {
    # maximum (t1) and minimum (t2) colonization time: between the present
    # and t2 the colonization time is fixed
    system2 <- DAISIE_DE_DD_system("interval2_ES", 1, rates, n_cap)
    y <- DAISIE_DE_DD_state(system2, initial)
    y <- DAISIE_DE_DD_solve_branch(system2, y, c(0, t2), methode, abstolint, reltolint)
    initial <- list(DA3 = DAISIE_DE_DD_get(system2, y, "DA3"),
                    DM3 = DAISIE_DE_DD_get(system2, y, "DM3"),
                    DM2 = DAISIE_DE_DD_get(system2, y, "DM2"))
    t_start <- t2
  }
  logL <- DAISIE_DE_DD_stem_to_island_age(DA3 = initial$DA3, DM3 = initial$DM3,
                                          DM2 = initial$DM2,
                                          t_start = t_start, t1 = t1, t0 = t0,
                                          stac = stac_stem, rates = rates,
                                          n_cap = n_cap, methode = methode,
                                          atol = abstolint, rtol = reltolint)
  logL + DAISIE_DE_DD_sampling_term(missnumspec, S, sampling)
}

DAISIE_DE_DD_logpEC <- function(brts,
                                missnumspec = 0,
                                pars1,
                                pars2 = c(100, 11, 0, 0),
                                stac,
                                methode = "odeint::runge_kutta_cash_karp54",
                                reltolint = 1e-10,
                                abstolint = 1e-16,
                                sampling = 'n',
                                dd_recolonization = FALSE,
                                mu_NE = pars1[2]) {
  if (!(stac %in% c(2, 3, 6, 7))) {
    stop("stac must be 2, 3, 6 or 7 for this function.")
  }
  t0 <- brts[1]
  t1 <- brts[2]
  # branching times from the youngest to the crown (t2)
  ti <- sort(brts)
  ti <- ti[1:(length(ti) - 2)]
  S <- length(brts) - 1
  n_cap <- pars2[1] + missnumspec
  rates <- DAISIE_DE_DD_rates(pars1, pars2[2], dd_recolonization, mu_NE)
  hidden <- DAISIE_DE_DD_hidden_init(missnumspec, S, sampling, n_cap)
  if (stac %in% c(3, 7)) {
    values <- list(DM3 = hidden)
  } else {
    values <- list(DA3 = hidden)
  }
  # k lineages of the observed clade; going back in time, k drops by one at
  # each branching time
  k <- S
  t_prev <- 0
  for (idx in seq_along(ti)) {
    system1 <- DAISIE_DE_DD_system("interval2_EC", k, rates, n_cap)
    y <- DAISIE_DE_DD_state(system1, values)
    y <- DAISIE_DE_DD_solve_branch(system1, y, c(t_prev, ti[idx]), methode,
                                   abstolint, reltolint)
    DA3 <- DAISIE_DE_DD_get(system1, y, "DA3")
    DM3 <- DAISIE_DE_DD_get(system1, y, "DM3")
    n_A <- seq_along(DA3) - 1
    n_M <- seq_along(DM3) - 1
    # speciation of one of the k - 1 lineages present just before the
    # branching time, at the rate for the number of species at that moment
    values <- list(DA3 = rates$clado((k - 1) + n_A) * DA3,
                   DM3 = rates$clado((k - 1) + n_M + 1) * DM3)
    if (idx == length(ti)) {
      # crown: it can also be the cladogenesis of the non-endemic colonist
      # (k = 0 lineages before it), which leaves the two observed lineages
      values$DM2 <- rates$clado(n_A + 1) * DA3
    }
    k <- k - 1
    t_prev <- ti[idx]
  }
  stac_stem <- if (stac %in% c(2, 3)) 2 else 6
  logL <- DAISIE_DE_DD_stem_to_island_age(DA3 = values$DA3, DM3 = values$DM3,
                                          DM2 = values$DM2,
                                          t_start = t_prev, t1 = t1, t0 = t0,
                                          stac = stac_stem, rates = rates,
                                          n_cap = n_cap, methode = methode,
                                          atol = abstolint, rtol = reltolint)
  logL + DAISIE_DE_DD_sampling_term(missnumspec, S, sampling)
}
