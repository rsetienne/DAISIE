# Diversity-dependent version of the DAISIE_DE equations.
#
# In DAISIE_DE, E describes the hidden species (species that are not in the
# data) one lineage at a time. With diversity dependence the rates depend on
# the total number of species N of the clade on the island, so the lineages
# are no longer independent. Instead of E, every DAISIE_DE variable (DA1,
# DM1, DM2, DM3, DA2, DA3) now carries an index n, the number of hidden
# species alive at that time. DA3 and DM3 also contain the k lineages of the
# observed clade (in DAISIE_DE these were the factor DE), so
#   DA3[n] (k lineages) corresponds to DE^k * DA3 * E^n when K = Inf,
#   DM3[n] (k lineages) corresponds to DE^k * DM3 * E^n,
#   DM2[n], DM1[n], DA2[n], DA1[n] correspond to DM2 * E^n, etc.
# The number of species is N = k + n (+ 1 if the mainland species is present
# on the island as a non-endemic species).
#
# The equations are linear in the variables, dy/dt = Q y. The functions in
# this file build the matrix Q for each interval of DAISIE_DE (same interval
# names as in DAISIE_DE_intervals.R) and integrate the system.
#
# pars1 (as in DAISIE_loglik_CS):
# 1: lambda_c
# 2: mu (endemic species; the non-endemic species uses mu_NE)
# 3: K
# 4: gamma
# 5: lambda_a

# Rates per species when the clade has N species on the island.
# ddep as in DAISIE: 0 no diversity dependence, 1 linear in cladogenesis,
# 11 linear in cladogenesis and colonization, 2 exponential in cladogenesis,
# 21 exponential in cladogenesis and colonization.
DAISIE_DE_DD_rates <- function(pars1,
                               ddep = 11,
                               dd_recolonization = FALSE,
                               mu_NE = pars1[2]) {
  if (!(ddep %in% c(0, 1, 11, 2, 21))) {
    stop("DAISIE_DE_DD: ddep must be 0, 1, 11, 2 or 21.")
  }
  lac <- pars1[1]
  K   <- pars1[3]
  gam <- pars1[4]
  clado <- function(N) {
    if (ddep == 1 || ddep == 11) return(pmax(0, lac * (1 - N / K)))
    if (ddep == 2 || ddep == 21) return(lac * exp(-N / K))
    return(rep(lac, length(N)))
  }
  colon <- function(N) {
    if (ddep == 11) return(pmax(0, gam * (1 - N / K)))
    if (ddep == 21) return(gam * exp(-N / K))
    return(rep(gam, length(N)))
  }
  # recolonization: colonization by the mainland species while it is already
  # present on the island as a non-endemic species (resets the colonization
  # time, adds no species); as in DAISIE_loglik_CS (dd_recolonization)
  recol <- function(N) {
    if (dd_recolonization) return(colon(N))
    return(rep(gam, length(N)))
  }
  list(clado = clado,
       colon = colon,
       recol = recol,
       mu_E = pars1[2],
       mu_NE = mu_NE,
       laa = pars1[5],
       K = K,
       ddep = ddep)
}

# Largest number of hidden species kept for a variable with k lineages of the
# observed clade and f = 1 if the non-endemic species is present (else f = 0).
# With linear diversity dependence in both cladogenesis and colonization
# (ddep = 11) the clade can never have more than ceiling(K) species, so this
# limit is exact. Otherwise the list is cut at n_cap (as lx in DAISIE).
DAISIE_DE_DD_n_max <- function(k, f, rates, n_cap) {
  if (rates$ddep == 11 && is.finite(rates$K)) {
    return(min(ceiling(rates$K) - k - f, n_cap))
  }
  return(n_cap)
}

# The variables of each interval: name = c(k, f)
DAISIE_DE_DD_blocks <- function(interval, k) {
  switch(interval,
         interval4    = list(DA1 = c(0, 0), DM1 = c(0, 1)),
         interval2_NE = list(DM2 = c(0, 1)),
         interval3_NE = list(DM1 = c(0, 1), DM2 = c(0, 1), DA2 = c(0, 0)),
         interval2_ES = list(DA3 = c(k, 0), DM3 = c(k, 1), DM2 = c(0, 1)),
         interval3_ES = list(DA3 = c(k, 0), DM3 = c(k, 1), DM1 = c(0, 1),
                             DM2 = c(0, 1), DA2 = c(0, 0)),
         interval2_EC = list(DA3 = c(k, 0), DM3 = c(k, 1)),
         stop("Unknown interval: ", interval))
}

# Build the system dy/dt = Q y of one interval.
# k: number of lineages of the observed clade in this interval.
DAISIE_DE_DD_system <- function(interval, k, rates, n_cap) {
  blocks <- DAISIE_DE_DD_blocks(interval, k)
  # position of each variable in the state vector
  layout <- list()
  offset <- 0
  for (b in names(blocks)) {
    kb <- blocks[[b]][1]
    fb <- blocks[[b]][2]
    n_max <- DAISIE_DE_DD_n_max(kb, fb, rates, n_cap)
    n <- if (n_max >= 0) 0:n_max else integer(0)
    layout[[b]] <- list(k = kb, f = fb, n = n, N = kb + n + fb,
                        index = offset + seq_along(n))
    offset <- offset + length(n)
  }
  rows <- c(); cols <- c(); vals <- c()

  # loss: the variable leaves its state at this rate
  loss <- function(b, rate) {
    v <- layout[[b]]
    rate <- rep_len(rate, length(v$n))
    rows <<- c(rows, v$index); cols <<- c(cols, v$index); vals <<- c(vals, -rate)
  }
  # gain: the variable b at n receives rate * (variable to at n + shift)
  gain <- function(b, to, rate, shift = 0) {
    v <- layout[[b]]; w <- layout[[to]]
    rate <- rep_len(rate, length(v$n))
    target <- v$n + shift
    ok <- target >= 0 & target <= max(c(-1, w$n)) & rate != 0
    rows <<- c(rows, v$index[ok])
    cols <<- c(cols, w$index[match(target[ok], w$n)])
    vals <<- c(vals, rate[ok])
  }

  for (b in names(layout)) {
    v <- layout[[b]]
    if (length(v$n) == 0) next
    lam <- rates$clado(v$N)
    # hidden species: speciation (n -> n + 1) and extinction (n -> n - 1)
    loss(b, (lam + rates$mu_E) * v$n)
    gain(b, b, lam * v$n, +1)
    gain(b, b, rates$mu_E * v$n, -1)
    # lineages of the observed clade: speciation with a hidden daughter
    # (either daughter can be the hidden one, hence 2); extinction not allowed
    loss(b, (lam + rates$mu_E) * v$k)
    gain(b, b, 2 * lam * v$k, +1)
  }

  # the non-endemic species disappears from the island and leaves 0, 1 or 2
  # hidden endemic species: extinction, anagenesis, cladogenesis
  ne_to <- function(b, to) {
    v <- layout[[b]]
    lam <- rates$clado(v$N)
    loss(b, rates$mu_NE + rates$laa + lam)
    gain(b, to, rates$mu_NE, 0)
    gain(b, to, rates$laa, +1)
    gain(b, to, lam, +2)
  }
  # the absent mainland species colonizes the island
  colonize <- function(b, to) {
    v <- layout[[b]]
    g <- rates$colon(v$N)
    loss(b, g)
    gain(b, to, g, 0)
  }

  if (interval == "interval4") {
    colonize("DA1", "DM1")
    ne_to("DM1", "DA1")
  } else if (interval == "interval2_NE") {
    # the observed non-endemic species stays as it is until the present:
    # every change, and recolonization (which would reset the colonization
    # time), is not allowed
    v <- layout$DM2
    loss("DM2", rates$clado(v$N) + rates$mu_NE + rates$laa + rates$recol(v$N))
  } else if (interval == "interval3_NE") {
    ne_to("DM1", "DA2")
    v <- layout$DM1
    loss("DM1", rates$recol(v$N))
    gain("DM1", "DM2", rates$recol(v$N), 0)
    ne_to("DM2", "DA2")
    colonize("DA2", "DM2")
  } else if (interval %in% c("interval2_ES", "interval3_ES", "interval2_EC")) {
    colonize("DA3", "DM3")
    ne_to("DM3", "DA3")
  }
  if (interval %in% c("interval2_ES", "interval3_ES")) {
    # DM2: the non-endemic colonist that becomes the observed lineage
    v <- layout$DM2
    lam <- rates$clado(v$N)
    loss("DM2", rates$laa + lam + rates$mu_NE)
    gain("DM2", "DA3", rates$laa, 0)        # anagenesis: becomes the observed lineage
    gain("DM2", "DA3", 2 * lam, +1)         # cladogenesis: observed + hidden species
    if (interval == "interval2_ES") {
      # the colonization time is known: recolonization is not allowed
      loss("DM2", rates$recol(v$N))
    } else {
      # maximum colonization time: the colonist may also leave hidden
      # species only, and the observed colonization comes later
      gain("DM2", "DA2", rates$mu_NE, 0)
      gain("DM2", "DA2", rates$laa, +1)
      gain("DM2", "DA2", lam, +2)
      ne_to("DM1", "DA2")
      w <- layout$DM1
      loss("DM1", rates$recol(w$N))
      gain("DM1", "DM2", rates$recol(w$N), 0)
      colonize("DA2", "DM2")
    }
  }
  list(layout = layout,
       size = offset,
       i = rows - 1L, j = cols - 1L, x = vals)
}

# Empty state vector of a system, filled from a named list of variables.
# Each element of values is a vector over n = 0, 1, ...; values that do not
# fit (n beyond the limit) are dropped, missing values are 0.
DAISIE_DE_DD_state <- function(system, values = list()) {
  y <- rep(0, system$size)
  for (b in names(values)) {
    v <- system$layout[[b]]
    m <- min(length(v$n), length(values[[b]]))
    if (m > 0) y[v$index[1:m]] <- values[[b]][1:m]
  }
  y
}

# One variable (vector over n) from a state vector
DAISIE_DE_DD_get <- function(system, y, b) {
  y[system$layout[[b]]$index]
}

# Integrate dy/dt = Q y over time = c(start, end); returns the state at the end
DAISIE_DE_DD_solve_branch <- function(system,
                                      initial_conditions,
                                      time,
                                      methode = "odeint::runge_kutta_cash_karp54",
                                      atol = 1e-16,
                                      rtol = 1e-10) {
  if (time[1] == time[2] || system$size == 0) return(initial_conditions)
  if (startsWith(methode, "odeint::")) {
    return(.Call("DAISIE_DE_DD_cpp_solve",
                 as.integer(system$i), as.integer(system$j),
                 as.numeric(system$x),
                 as.numeric(initial_conditions), as.numeric(time),
                 methode, atol, rtol))
  }
  Q <- Matrix::sparseMatrix(i = system$i + 1, j = system$j + 1,
                            x = system$x, dims = c(system$size, system$size))
  rhs <- function(t, y, parms) list(as.vector(Q %*% y))
  solution <- deSolve::ode(y = initial_conditions,
                           times = time,
                           func = rhs,
                           parms = NULL,
                           method = methode,
                           atol = atol,
                           rtol = rtol)
  as.vector(solution[2, -1])
}
