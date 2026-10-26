#' Expected numbers of endemic and non-endemic species per range in the
#' archipelago model without diversity dependence
#'
#' @description
#' Generalisation of \code{\link{DAISIE_ExpEIN}()} to the archipelago model
#' DAISIE_A. For infinite carrying capacities lineages do not interact, and
#' the expected numbers of species in every range state satisfy a linear
#' system of ordinary differential equations:
#' \itemize{
#'   \item for every mainland species, the probability \eqn{p_R(t)} that it is
#'     present as a non-endemic lineage with range \eqn{R} (and \eqn{p_0(t)}
#'     that it is not present as a non-endemic) follows the forward Kolmogorov
#'     equation of a continuous-time Markov chain on
#'     \eqn{\{0\} \cup \mathcal{R}};
#'   \item the expected number \eqn{e_R(t)} of endemic species with range
#'     \eqn{R} descending from that mainland species follows the mean
#'     equation of a multi-type branching process with immigration from the
#'     non-endemic chain.
#' }
#' The expected numbers for the whole mainland pool are \eqn{M p_R(t)} and
#' \eqn{M e_R(t)}.
#'
#' With \code{lumped = TRUE} the lumped (aggregated) system for exchangeable
#' islands is solved, in which the state of a lineage is only the number of
#' islands it occupies. This reduces the dimension from \eqn{2^n - 1} to
#' \eqn{n} (or \code{max_range_size}) range states, and is exact when all
#' islands have the same rates, all pairs of islands are equally connected,
#' and all islands are equally old.
#'
#' @param tvec Vector of times (since the birth of the oldest island) at which
#'   the expectations are computed.
#' @param pars A \code{DAISIE_A_pars} object (all \code{K} are ignored and
#'   treated as infinite).
#' @param M Number of mainland species.
#' @param lumped Logical; solve the lumped system for exchangeable islands.
#' @param total_time Age of the oldest island; only needed when
#'   \code{pars$island_ages} is not \code{NULL} (default \code{max(tvec)}).
#' @param rtol,atol Tolerances of the ODE solver.
#'
#' @return A list with the times \code{tvec}, the expected numbers of endemic,
#'   non-endemic and all species on the archipelago (\code{ExpE},
#'   \code{ExpI}, \code{ExpN}), per island (\code{ExpE_island},
#'   \code{ExpI_island}, \code{ExpN_island}, \code{ExpSIE_island} for
#'   single-island endemics; matrices with one column per island), and per
#'   range state (\code{ExpE_range}, \code{ExpI_range}; one column per range,
#'   named by the range string, or by the range size if \code{lumped}).
#' @author Rampal S. Etienne
#' @family DAISIE_A
#' @export
#' @examples
#' pars <- DAISIE_A_create_pars(n_islands = 3, lac = 0.5, mu = 0.2, gam = 0.02,
#'   laa = 0.5, lav = 0.2, disp = 0.3, connectivity = "line")
#' exp_vals <- DAISIE_A_ExpEIN(tvec = c(1, 2, 5), pars = pars, M = 100)
#' exp_vals$ExpN_island
DAISIE_A_ExpEIN <- function(tvec,
                            pars,
                            M,
                            lumped = FALSE,
                            total_time = max(tvec),
                            rtol = 1e-10,
                            atol = 1e-12) {
  if (!DAISIE_A_are_pars(pars)) {
    stop("'pars' must be created with DAISIE_A_create_pars().")
  }
  if (any(is.finite(pars$K))) {
    warning("K must be infinite; proceeding with infinite K.")
  }
  if (lumped) {
    return(DAISIE_A_ExpEIN_lumped(tvec = tvec, pars = pars, M = M,
                                  rtol = rtol, atol = atol))
  }
  n <- pars$n_islands
  ranges <- DAISIE_A_ranges(pars)
  S <- length(ranges$sets)
  emerge <- DAISIE_A_emergence_times(pars, total_time)
  submerge <- DAISIE_A_submergence_times(pars, total_time)
  breaks <- sort(unique(c(emerge[emerge > 0], submerge[is.finite(submerge)])))
  x0 <- c(1, rep(0, 2 * S))
  tout <- sort(unique(c(0, tvec)))
  bounds <- c(0, breaks[breaks < max(tout)], max(tout))
  res <- matrix(NA_real_, nrow = length(tout), ncol = length(x0))
  res[1, ] <- x0
  x <- x0
  for (b in seq_len(length(bounds) - 1)) {
    # islands submerging now lose all their populations
    for (i in which(submerge == bounds[b])) {
      x <- DAISIE_A_submergence_map(x, i, ranges)
    }
    exists <- as.numeric(emerge <= bounds[b] & bounds[b] < submerge)
    G <- DAISIE_A_generator(pars = pars, ranges = ranges, exists = exists)
    tt <- tout[tout > bounds[b] & tout <= bounds[b + 1]]
    times <- unique(c(bounds[b], tt, bounds[b + 1]))
    if (length(times) < 2) next
    sol <- deSolve::ode(
      y = x,
      times = times,
      func = function(t, y, parms) list(as.vector(parms %*% y)),
      parms = G,
      method = "lsoda",
      rtol = rtol,
      atol = atol
    )
    sol <- sol[, -1, drop = FALSE]
    for (t in tt) {
      res[which(tout == t), ] <- sol[which(times == t), ]
    }
    x <- sol[nrow(sol), ]
  }
  res <- res[match(tvec, tout), , drop = FALSE]
  p <- res[, 1 + seq_len(S), drop = FALSE]
  e <- res[, 1 + S + seq_len(S), drop = FALSE]
  colnames(p) <- colnames(e) <- ranges$strings
  member <- matrix(0, nrow = S, ncol = n)
  for (s in seq_len(S)) member[s, ranges$sets[[s]]] <- 1
  single <- which(vapply(ranges$sets, length, numeric(1)) == 1)
  sie <- matrix(0, nrow = length(tvec), ncol = n)
  for (s in single) sie[, ranges$sets[[s]]] <- e[, s]
  out <- list(
    tvec = tvec,
    ExpE = M * rowSums(e),
    ExpI = M * rowSums(p),
    ExpN = M * (rowSums(e) + rowSums(p)),
    ExpE_island = M * e %*% member,
    ExpI_island = M * p %*% member,
    ExpN_island = M * (e + p) %*% member,
    ExpSIE_island = M * sie,
    ExpE_range = M * e,
    ExpI_range = M * p
  )
  for (nm in c("ExpE_island", "ExpI_island", "ExpN_island", "ExpSIE_island")) {
    colnames(out[[nm]]) <- pars$island_names
  }
  return(out)
}

#' Effect of the submergence of an island on the expected numbers
#'
#' @param x State vector \eqn{(p_0, p_R, e_R)} of the expected-value system.
#' @param i Index of the submerging island.
#' @param ranges Output of \code{\link{DAISIE_A_ranges}()}.
#' @return The state vector after island \eqn{i} has lost all its
#'   populations: \eqn{p_R} and \eqn{e_R} move to \eqn{R \setminus \{i\}},
#'   \eqn{p_{\{i\}}} moves to \eqn{p_0} and \eqn{e_{\{i\}}} is lost.
#' @keywords internal
DAISIE_A_submergence_map <- function(x, i, ranges) {
  S <- length(ranges$sets)
  y <- x
  hit <- which(vapply(ranges$sets, function(R) i %in% R, logical(1)))
  y[1 + hit] <- 0
  y[1 + S + hit] <- 0
  for (s in hit) {
    R <- ranges$sets[[s]]
    if (length(R) == 1) {
      y[1] <- y[1] + x[1 + s]
    } else {
      t <- match(sum(2^(setdiff(R, i) - 1)), ranges$codes)
      if (is.na(t)) stop("Range after submergence is not an allowed range.")
      y[1 + t] <- y[1 + t] + x[1 + s]
      y[1 + S + t] <- y[1 + S + t] + x[1 + S + s]
    }
  }
  return(y)
}

#' Generator of the linear system for the expected numbers of species
#'
#' @param pars A \code{DAISIE_A_pars} object.
#' @param ranges Output of \code{\link{DAISIE_A_ranges}()}.
#' @param exists 0/1 vector indicating which islands have emerged.
#'
#' @return Matrix \eqn{G} such that \eqn{dx/dt = G x} with
#'   \eqn{x = (p_0, p_R, e_R)}.
#' @keywords internal
DAISIE_A_generator <- function(pars, ranges, exists) {
  n <- pars$n_islands
  S <- length(ranges$sets)
  code_index <- integer(2^n - 1)
  code_index[ranges$codes] <- seq_len(S)
  idx <- function(set) code_index[sum(2^(set - 1))]
  gam <- pars$gam * exists
  D <- sweep(pars$D, 2, exists, "*")
  laa_isl <- DAISIE_A_laa_isl(pars)
  Q <- matrix(0, S + 1, S + 1)
  B <- matrix(0, S, S)
  A <- matrix(0, S, S)
  for (i in seq_len(n)) {
    Q[1, 1 + idx(i)] <- Q[1, 1 + idx(i)] + gam[i]
  }
  for (s in seq_len(S)) {
    R <- ranges$sets[[s]]
    r <- length(R)
    out_isl <- setdiff(seq_len(n), R)
    # range expansion: recolonisation (non-endemics only) and dispersal
    for (j in out_isl) {
      Rj <- sort(c(R, j))
      if (DAISIE_A_range_allowed(Rj, pars)) {
        sj <- idx(Rj)
        Q[1 + s, 1 + sj] <- Q[1 + s, 1 + sj] + gam[j] + sum(D[R, j])
        A[s, s] <- A[s, s] - sum(D[R, j])
        A[s, sj] <- A[s, sj] + sum(D[R, j])
      }
    }
    # local and global extinction
    if (r == 1) {
      Q[1 + s, 1] <- Q[1 + s, 1] + pars$mu[R]
      A[s, s] <- A[s, s] - pars$mu[R]
    } else {
      for (i in R) {
        Ri <- setdiff(R, i)
        if (DAISIE_A_range_allowed(Ri, pars)) {
          si <- idx(Ri)
          Q[1 + s, 1 + si] <- Q[1 + s, 1 + si] + pars$mu[i]
          A[s, s] <- A[s, s] - pars$mu[i]
          A[s, si] <- A[s, si] + pars$mu[i]
        }
      }
    }
    # anagenesis
    Q[1 + s, 1] <- Q[1 + s, 1] + pars$laa
    B[s, s] <- B[s, s] + pars$laa
    # in-situ cladogenesis (split or budding)
    for (i in R) {
      si <- idx(i)
      if (r == 1 && pars$clado_single == "split") {
        Q[1 + s, 1] <- Q[1 + s, 1] + pars$lac[i]
        B[s, si] <- B[s, si] + 2 * pars$lac[i]
      } else {
        B[s, si] <- B[s, si] + pars$lac[i]
      }
      A[s, si] <- A[s, si] + pars$lac[i]
    }
    # vicariance
    if (r >= 2 && pars$lav > 0) {
      # rate of every admissible split and total vicariance rate of R
      vr <- DAISIE_A_vic_rates(R, pars)
      Q[1 + s, 1] <- Q[1 + s, 1] + vr$total
      A[s, s] <- A[s, s] - vr$total
      for (k in seq_along(vr$splits)) {
        sa <- idx(sort(vr$splits[[k]]))
        sb <- idx(sort(setdiff(R, vr$splits[[k]])))
        B[s, sa] <- B[s, sa] + vr$rates[k]
        B[s, sb] <- B[s, sb] + vr$rates[k]
        A[s, sa] <- A[s, sa] + vr$rates[k]
        A[s, sb] <- A[s, sb] + vr$rates[k]
      }
    }
    # island-level anagenesis of non-endemics: the population on island i
    # becomes a new endemic on {i}, the non-endemic keeps R \ {i}
    if (r >= 2) {
      for (i in R) {
        Ri <- setdiff(R, i)
        if (laa_isl[i] > 0 && DAISIE_A_range_allowed(Ri, pars)) {
          si <- idx(Ri)
          Q[1 + s, 1 + si] <- Q[1 + s, 1 + si] + laa_isl[i]
          B[s, idx(i)] <- B[s, idx(i)] + laa_isl[i]
        }
      }
    }
  }
  diag(Q) <- 0
  diag(Q) <- -rowSums(Q)
  G <- matrix(0, 1 + 2 * S, 1 + 2 * S)
  G[1:(S + 1), 1:(S + 1)] <- t(Q)
  G[S + 1 + seq_len(S), 1 + seq_len(S)] <- t(B)
  G[S + 1 + seq_len(S), S + 1 + seq_len(S)] <- t(A)
  return(G)
}

#' Expected numbers of species for exchangeable islands (lumped system)
#'
#' @inheritParams DAISIE_A_ExpEIN
#' @return See \code{\link{DAISIE_A_ExpEIN}()}.
#' @keywords internal
DAISIE_A_ExpEIN_lumped <- function(tvec, pars, M, rtol = 1e-10,
                                   atol = 1e-12) {
  n <- pars$n_islands
  C <- pars$C
  offdiag <- C[row(C) != col(C)]
  same <- function(x) all(abs(x - x[1]) <= 1e-12 * max(1, abs(x[1])))
  laa_isl <- DAISIE_A_laa_isl(pars)
  if (!(same(pars$lac) && same(pars$mu) && same(pars$gam) &&
        same(laa_isl) &&
        (n == 1 || same(offdiag)) &&
        (is.null(pars$island_ages) || same(pars$island_ages)))) {
    stop("The lumped system is only exact for exchangeable islands: equal ",
         "lac, mu, gam, laa_isl and island ages on all islands and equal ",
         "connectivity between all pairs of islands.")
  }
  if (!is.null(pars$submergence_ages)) {
    stop("The lumped system does not support island submergence.")
  }
  if (!is.null(pars$island_ages) && max(tvec) > pars$island_ages[1]) {
    stop("All times must be within the age of the islands.")
  }
  lac <- pars$lac[1]
  mu <- pars$mu[1]
  gam <- pars$gam[1]
  d <- if (n > 1) pars$D[1, 2] else 0
  lia <- laa_isl[1]
  kmax <- pars$max_range_size
  # vicariance of a range of size k: expected numbers of daughters of each
  # size per event and total weight W (all pairs of islands have
  # connectivity cc, so a split into parts of sizes j and k - j cuts
  # 2 * cc * j * (k - j))
  cc <- if (n > 1) C[1, 2] else 0
  kappa <- DAISIE_A_vic_kappa(pars)
  vic_lumped <- function(k) {
    out <- rep(0, kmax)
    if (pars$vic_mode == "single") {
      if (k == 2) {
        out[1] <- 2
      } else {
        out[1] <- 1
        out[k - 1] <- out[k - 1] + 1
      }
      W <- (if (k == 2) 1 else k) * exp(-kappa * 2 * cc * (k - 1))
    } else {
      j <- seq_len(k - 1)
      wj <- exp(-kappa * 2 * cc * j * (k - j))
      cnt <- choose(k, j)
      W <- sum(cnt * wj) / 2
      out[j] <- out[j] + cnt * wj / W
    }
    list(daughters = out, W = W)
  }
  Q <- matrix(0, kmax + 1, kmax + 1)
  B <- matrix(0, kmax, kmax)
  A <- matrix(0, kmax, kmax)
  Q[1, 2] <- n * gam
  for (k in seq_len(kmax)) {
    up <- k * (n - k) * d
    if (k < kmax) {
      Q[1 + k, 2 + k] <- (n - k) * gam + up
      A[k, k] <- A[k, k] - up
      A[k, k + 1] <- A[k, k + 1] + up
    }
    if (k == 1) {
      Q[2, 1] <- Q[2, 1] + mu
      A[1, 1] <- A[1, 1] - mu
    } else {
      Q[1 + k, k] <- Q[1 + k, k] + k * mu
      A[k, k] <- A[k, k] - k * mu
      A[k, k - 1] <- A[k, k - 1] + k * mu
    }
    Q[1 + k, 1] <- Q[1 + k, 1] + pars$laa
    B[k, k] <- B[k, k] + pars$laa
    if (k == 1 && pars$clado_single == "split") {
      Q[2, 1] <- Q[2, 1] + lac
      B[1, 1] <- B[1, 1] + 2 * lac
    } else {
      B[k, 1] <- B[k, 1] + k * lac
    }
    A[k, 1] <- A[k, 1] + k * lac
    if (k >= 2 && pars$lav > 0) {
      vl <- vic_lumped(k)
      vd <- vl$daughters
      vic_total <- if (DAISIE_A_vic_per(pars) == "split") {
        pars$lav * vl$W
      } else {
        pars$lav
      }
      Q[1 + k, 1] <- Q[1 + k, 1] + vic_total
      A[k, k] <- A[k, k] - vic_total
      B[k, ] <- B[k, ] + vic_total * vd
      A[k, ] <- A[k, ] + vic_total * vd
    }
    if (k >= 2) {
      # island-level anagenesis of non-endemics on one of the k islands
      ia <- k * lia
      Q[1 + k, k] <- Q[1 + k, k] + ia
      B[k, 1] <- B[k, 1] + ia
    }
  }
  diag(Q) <- 0
  diag(Q) <- -rowSums(Q)
  G <- matrix(0, 1 + 2 * kmax, 1 + 2 * kmax)
  G[1:(kmax + 1), 1:(kmax + 1)] <- t(Q)
  G[kmax + 1 + seq_len(kmax), 1 + seq_len(kmax)] <- t(B)
  G[kmax + 1 + seq_len(kmax), kmax + 1 + seq_len(kmax)] <- t(A)
  tout <- sort(unique(c(0, tvec)))
  sol <- deSolve::ode(
    y = c(1, rep(0, 2 * kmax)),
    times = tout,
    func = function(t, y, parms) list(as.vector(parms %*% y)),
    parms = G,
    method = "lsoda",
    rtol = rtol,
    atol = atol
  )
  res <- sol[match(tvec, tout), -1, drop = FALSE]
  p <- res[, 1 + seq_len(kmax), drop = FALSE]
  e <- res[, 1 + kmax + seq_len(kmax), drop = FALSE]
  colnames(p) <- colnames(e) <- paste0("size_", seq_len(kmax))
  frac <- seq_len(kmax) / n
  isl <- function(x) matrix(rep(x, n), ncol = n,
                            dimnames = list(NULL, pars$island_names))
  out <- list(
    tvec = tvec,
    ExpE = M * rowSums(e),
    ExpI = M * rowSums(p),
    ExpN = M * (rowSums(e) + rowSums(p)),
    ExpE_island = isl(M * as.vector(e %*% frac)),
    ExpI_island = isl(M * as.vector(p %*% frac)),
    ExpN_island = isl(M * as.vector((e + p) %*% frac)),
    ExpSIE_island = isl(M * e[, 1] / n),
    ExpE_range = M * e,
    ExpI_range = M * p
  )
  return(out)
}
