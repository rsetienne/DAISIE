#' @title Truncated power series for n-sampling in DAISIE_DE
#' @description Internal helpers for the power-series method of
#'   \code{DAISIE_DE_n}. Under n-sampling the likelihood of a clade with
#'   \code{S} sampled and \code{m} missing species is
#'   \code{[z^m] L(z) / choose(S + m, S)}, where \code{L(z)} is the
#'   likelihood of the D-E equations with initial condition \code{E(0) = z}
#'   (and \code{D(0) = 1} for every sampled tip). Every state variable is
#'   therefore a power series in \code{z}. Because the D-E equations are
#'   polynomial in the state variables, the coefficients of these series
#'   satisfy a closed system of ordinary differential equations (obtained by
#'   Cauchy products), in which the coefficient of degree \code{k} only
#'   depends on coefficients of degree \code{<= k}. Truncating after degree
#'   \code{m} is therefore exact, and the likelihood is obtained from a single
#'   solve of a system that is \code{m + 1} times larger than the original
#'   one, with real arithmetic only.
#' @name DAISIE_DE_series
#' @keywords internal
#' @noRd
NULL

#' Truncated Cauchy product of two coefficient vectors
#' @noRd
DAISIE_DE_series_conv <- function(a, b) {
  P <- length(a)
  out <- numeric(P)
  for (k in seq_len(P)) {
    out[k] <- sum(a[seq_len(k)] * b[k:1])
  }
  return(out)
}

#' Builds a state matrix (rows: coefficients of degree 0..P-1, columns:
#' variables) from scalars (constants) and coefficient vectors.
#' @noRd
DAISIE_DE_series_state <- function(P, ...) {
  vals <- list(...)
  out <- matrix(0, nrow = P, ncol = length(vals),
                dimnames = list(NULL, names(vals)))
  for (v in names(vals)) {
    x <- vals[[v]]
    if (length(x) == 1) {
      out[1, v] <- x
    } else {
      out[, v] <- x
    }
  }
  return(out)
}

#' Coefficient vector of the identity z, truncated after degree P - 1
#' @noRd
DAISIE_DE_series_z <- function(P) {
  z <- numeric(P)
  if (P >= 2) z[2] <- 1
  return(z)
}

#' Right-hand side of the coefficient equations for deSolve. The state
#' vector holds the coefficients variable by variable, in the same order as
#' in the corresponding interval functions.
#' @noRd
DAISIE_DE_series_rhs <- function(t, y, parms) {
  P <- parms$P
  lc <- parms$pars[1]
  mu_E <- parms$pars[2]
  mu_NE <- parms$pars[3]
  g <- parms$pars[4]
  la <- parms$pars[5]
  conv <- DAISIE_DE_series_conv
  x <- matrix(y, nrow = P)
  e1 <- c(1, rep(0, P - 1))
  dE <- function(E, EE) mu_E * e1 - (mu_E + lc) * E + lc * EE
  f_of <- function(E, EE) mu_NE * e1 + la * E + lc * EE
  dx <- switch(
    parms$interval,
    interval2_NE = {
      DM2 <- x[, 1]; E <- x[, 2]
      EE <- conv(E, E)
      cbind(-(lc + mu_NE + g + la) * DM2, dE(E, EE))
    },
    interval2_ES = {
      DE <- x[, 1]; DM2 <- x[, 2]; DM3 <- x[, 3]; E <- x[, 4]; DA3 <- x[, 5]
      EE <- conv(E, E)
      DEE <- conv(DE, E)
      f <- f_of(E, EE)
      h <- la * DE + 2 * lc * DEE
      cbind(-(lc + mu_E) * DE + 2 * lc * DEE,
            -(lc + mu_NE + g + la) * DM2 + conv(h, DA3),
            -(lc + mu_NE + la) * DM3 + conv(f, DA3),
            dE(E, EE),
            -g * DA3 + g * DM3)
    },
    interval2_EC = {
      DE <- x[, 1]; DM3 <- x[, 2]; E <- x[, 3]; DA3 <- x[, 4]
      EE <- conv(E, E)
      DEE <- conv(DE, E)
      f <- f_of(E, EE)
      cbind(-(lc + mu_E) * DE + 2 * lc * DEE,
            -(lc + mu_NE + la) * DM3 + conv(f, DA3),
            dE(E, EE),
            -g * DA3 + g * DM3)
    },
    interval3_ES = {
      DE <- x[, 1]; DM1 <- x[, 2]; DM2 <- x[, 3]; DM3 <- x[, 4]
      E <- x[, 5]; DA2 <- x[, 6]; DA3 <- x[, 7]
      EE <- conv(E, E)
      DEE <- conv(DE, E)
      f <- f_of(E, EE)
      h <- la * DE + 2 * lc * DEE
      fDA2 <- conv(f, DA2)
      cbind(-(lc + mu_E) * DE + 2 * lc * DEE,
            -(lc + mu_NE + la + g) * DM1 + g * DM2 + fDA2,
            -(lc + mu_NE + la) * DM2 + fDA2 + conv(h, DA3),
            -(lc + mu_NE + la) * DM3 + conv(f, DA3),
            dE(E, EE),
            -g * DA2 + g * DM2,
            -g * DA3 + g * DM3)
    },
    interval3_NE = {
      DM1 <- x[, 1]; DM2 <- x[, 2]; E <- x[, 3]; DA2 <- x[, 4]
      EE <- conv(E, E)
      fDA2 <- conv(f_of(E, EE), DA2)
      cbind(-(lc + mu_NE + la + g) * DM1 + fDA2 + g * DM2,
            -(lc + mu_NE + la) * DM2 + fDA2,
            dE(E, EE),
            -g * DA2 + g * DM2)
    },
    interval4 = {
      DA1 <- x[, 1]; DM1 <- x[, 2]; E <- x[, 3]
      EE <- conv(E, E)
      cbind(-g * DA1 + g * DM1,
            -(la + mu_NE + lc) * DM1 + conv(f_of(E, EE), DA1),
            dE(E, EE))
    },
    stop("Unknown interval: ", parms$interval)
  )
  return(list(as.vector(dx)))
}

#' Integrates the coefficient equations of one interval. Returns a list with
#' one state matrix (as produced by DAISIE_DE_series_state) per time point,
#' the first one being the initial state.
#' @noRd
DAISIE_DE_solve_branch_series <- function(interval_name,
                                          init,
                                          time,
                                          parameter,
                                          methode,
                                          atol,
                                          rtol) {
  P <- nrow(init)
  vars <- colnames(init)
  y0 <- as.vector(init)
  if (startsWith(methode, "odeint::")) {
    sol <- .Call("DAISIE_DE_cpp_solve_series",
                 parameter[[1]],     # lambda_c
                 parameter[[5]],     # lambda_a
                 parameter[[2]],     # mu_E
                 parameter[[3]],     # mu_NE
                 parameter[[4]],     # gamma
                 interval_name,
                 methode,
                 y0,
                 as.integer(P),
                 time,
                 atol,
                 rtol)
    states <- if (length(time) == 2) list(y0, sol$states) else sol$states
  } else {
    out <- deSolve::ode(y = y0,
                        times = time,
                        func = DAISIE_DE_series_rhs,
                        parms = list(pars = parameter,
                                     P = P,
                                     interval = interval_name),
                        method = methode,
                        atol = atol,
                        rtol = rtol)
    states <- lapply(seq_len(nrow(out)), function(i) out[i, -1])
  }
  return(lapply(states, function(s) {
    matrix(unlist(s), nrow = P, dimnames = list(NULL, vars))
  }))
}

#' Power-series version of DAISIE_DE_logpES: returns
#' log([z^missnumspec] L(z)), without the binomial normalisation.
#' @noRd
DAISIE_DE_logpES_series <- function(brts,
                                    missnumspec,
                                    pars1,
                                    stac,
                                    methode,
                                    reltolint,
                                    abstolint) {
  if (!(stac %in% c(2, 3, 5, 7, 9))) {
    stop("stac must be 2, 3, 5, 7 or 9 for this function.")
  }
  P <- missnumspec + 1
  t0 <- brts[1]
  t1 <- brts[2]
  t2 <- brts[3]
  tp <- 0
  z <- DAISIE_DE_series_z(P)
  st <- function(...) DAISIE_DE_series_state(P, ...)
  solve <- function(interval_name, init, time) {
    DAISIE_DE_solve_branch_series(interval_name = interval_name,
                                  init = init,
                                  time = time,
                                  parameter = pars1,
                                  methode = methode,
                                  atol = abstolint,
                                  rtol = reltolint)
  }

  interval_name <- "interval2_ES"
  time1 <- c(tp, t1)
  if (stac == 2) {
    init1 <- st(DE = 1, DM2 = 0, DM3 = 0, E = z, DA3 = 1)
  } else if (stac == 3) {
    init1 <- st(DE = 1, DM2 = 0, DM3 = 1, E = z, DA3 = 0)
  } else if (stac == 5) {
    init1 <- st(DE = 1, DM1 = 0, DM2 = 0, DM3 = 0, E = z, DA2 = 0, DA3 = 1)
    interval_name <- "interval3_ES"
  } else if (stac == 7) {
    init1 <- st(DE = 1, DM1 = 0, DM2 = 0, DM3 = 1, E = z, DA2 = 0, DA3 = 0)
    interval_name <- "interval3_ES"
  } else {
    # stac 9: as in DAISIE_DE_logpES, E = 0 at the present
    init1 <- st(DE = 1, DM2 = 0, DM3 = 0, E = 0, DA3 = 1)
    time1 <- c(tp, t2)
  }
  s1 <- solve(interval_name, init1, time1)[[2]]

  if (stac == 9) {
    init2 <- st(DE = s1[, "DE"], DM1 = 0, DM2 = s1[, "DM2"],
                DM3 = s1[, "DM3"], E = s1[, "E"], DA2 = 0, DA3 = s1[, "DA3"])
    s2 <- solve("interval3_ES", init2, c(t2, t1))[[2]]
    init3 <- st(DA1 = s2[, "DA2"], DM1 = s2[, "DM1"], E = s2[, "E"])
  } else if (stac == 5 || stac == 7) {
    init3 <- st(DA1 = s1[, "DA2"], DM1 = s1[, "DM1"], E = s1[, "E"])
  } else {
    init3 <- st(DA1 = pars1[4] * s1[, "DM2"],
                DM1 = pars1[4] * s1[, "DM2"],
                E = s1[, "E"])
  }
  s3 <- solve("interval4", init3, c(t1, t0))[[2]]
  return(log(unname(s3[P, "DA1"])))
}

#' Power-series version of DAISIE_DE_logpEC: returns
#' log([z^missnumspec] L(z)), without the binomial normalisation.
#' @noRd
DAISIE_DE_logpEC_series <- function(brts,
                                    missnumspec,
                                    stac,
                                    pars1,
                                    methode,
                                    reltolint,
                                    abstolint) {
  if (!(stac %in% c(2, 3, 6, 7))) {
    stop("stac must be 2, 3, 6 or 7 for this function.")
  }
  P <- missnumspec + 1
  t0 <- brts[1]
  t1 <- brts[2]
  t2 <- brts[3]
  ti <- sort(brts)
  ti <- ti[1:(length(ti) - 2)]
  z <- DAISIE_DE_series_z(P)
  conv <- DAISIE_DE_series_conv
  st <- function(...) DAISIE_DE_series_state(P, ...)
  solve <- function(interval_name, init, time) {
    DAISIE_DE_solve_branch_series(interval_name = interval_name,
                                  init = init,
                                  time = time,
                                  parameter = pars1,
                                  methode = methode,
                                  atol = abstolint,
                                  rtol = reltolint)
  }

  if (stac == 3 || stac == 7) {
    init1 <- st(DE = 1, DM3 = 1, E = z, DA3 = 0)
  } else {
    init1 <- st(DE = 1, DM3 = 0, E = z, DA3 = 1)
  }
  # a single tip branch evaluated at all branching times (comb)
  s0 <- solve("interval2_EC", init1, c(0, ti))
  times <- rbind(c(0, ti[1:(length(ti) - 1)]), ti)
  for (idx in seq_along(ti)) {
    s1 <- solve("interval2_EC", init1, times[, idx])[[2]]
    init1 <- st(DE = pars1[1] * conv(s0[[idx + 1]][, "DE"], s1[, "DE"]),
                DM3 = 0,
                E = s0[[idx + 1]][, "E"],
                DA3 = 1)
  }
  s0_end <- s0[[length(ti) + 1]]
  DE <- init1[, "DE"]
  DM2 <- conv(DE, s0_end[, "DA3"])
  if (stac == 6 || stac == 7) {
    init2 <- st(DE = DE, DM1 = 0, DM2 = DM2, DM3 = s0_end[, "DM3"],
                E = init1[, "E"], DA2 = 0, DA3 = s0_end[, "DA3"])
    s2 <- solve("interval3_ES", init2, c(t2, t1))[[2]]
    init3 <- st(DA1 = s2[, "DA2"], DM1 = s2[, "DM1"], E = s2[, "E"])
  } else {
    init2 <- st(DE = DE, DM2 = DM2, DM3 = s0_end[, "DM3"],
                E = init1[, "E"], DA3 = s0_end[, "DA3"])
    s2 <- solve("interval2_ES", init2, c(t2, t1))[[2]]
    init3 <- st(DA1 = pars1[4] * s2[, "DM2"],
                DM1 = pars1[4] * s2[, "DM2"],
                E = s2[, "E"])
  }
  s3 <- solve("interval4", init3, c(t1, t0))[[2]]
  return(log(unname(s3[P, "DA1"])))
}
