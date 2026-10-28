# Approximate likelihood of DAISIE_A for finite carrying capacities.
#
# Mean-field approximation: in the likelihood equations every rate of an
# event that adds a species to island i (colonisation of i, in-situ
# cladogenesis on i, dispersal into i) is multiplied by
#   f_i(t) = max(0, 1 - nu_i(t) / K_i),
# where nu_i(t) is a deterministic approximation of the number of species of
# the mainland species' descendants on island i at age t, given the data:
#   nu_i(t) = l_i(t) + h_i(t),
#   l_i(t)  lineages of the reconstructed phylogenies alive at age t, each
#           counted on island i with the fraction of its sampled descendants
#           whose range contains i (exact for one island);
#   h_i(t)  expected number of the other species (lineages without sampled
#           descendants: extinct before the present or missing from the
#           phylogeny) on island i, given the data, computed from the
#           branching process conditioned on leaving no sampled descendants.
# An unobserved lineage (equations for E) or the non-endemic lineage of a
# mainland species (equations for F) also counts itself when it speciates on
# its own island ("cavity" correction): its cladogenesis on island i is
# multiplied by f_i evaluated at nu_i + 1.
# Because the rates depend on nu and nu on the rates, nu is computed by
# fixed-point iteration. Given nu, the likelihood equations are those of the
# K = Inf model with time-dependent rates, solved as for K = Inf (including
# the power series for missing species). See the report for details.

#' Is diversity dependence active?
#'
#' @param pars A \code{DAISIE_A_pars} object.
#' @return \code{TRUE} if at least one carrying capacity is finite.
#' @keywords internal
DAISIE_A_dd_active <- function(pars) {
  any(is.finite(pars$K))
}

#' Split the rates into diversity-independent and per-island
#' diversity-dependent components
#'
#' @description The right-hand side of the likelihood equations is linear in
#'   the rates. Component 1 contains the rates that are not diversity
#'   dependent (extinction, anagenesis, island-level anagenesis, vicariance);
#'   component \eqn{1 + i} the cladogenesis rate on island \eqn{i} and
#'   component \eqn{1 + n + i} the rates of colonisation of and dispersal into
#'   island \eqn{i}; these are multiplied by diversity-dependence factors
#'   (see \code{\link{DAISIE_A_dd_wmat}()}).
#'
#' @param pars A \code{DAISIE_A_pars} object.
#' @return A list of \code{2 n_islands + 1} parameter objects.
#' @keywords internal
DAISIE_A_dd_components <- function(pars) {
  n <- pars$n_islands
  p0 <- pars
  p0$lac <- rep(0, n)
  p0$gam <- rep(0, n)
  p0$D <- 0 * pars$D
  base <- pars
  base$mu <- rep(0, n)
  base$laa <- 0
  base$laa_isl <- rep(0, n)
  base$lav <- 0
  base$lac <- rep(0, n)
  base$gam <- rep(0, n)
  base$D <- 0 * pars$D
  comps <- list(p0)
  for (i in seq_len(n)) {
    keep <- as.numeric(seq_len(n) == i)
    pc <- base
    pc$lac <- pars$lac * keep
    comps[[1 + i]] <- pc
    pd <- base
    pd$gam <- pars$gam * keep
    pd$D <- sweep(pars$D, 2, keep, "*")
    comps[[1 + n + i]] <- pd
  }
  comps
}

#' Weights of the rate components per block of equations
#'
#' @description With \eqn{\nu_i = l_i + h_i + p_i} (reconstructed lineages,
#'   unobserved species, presence of the non-endemic lineage of a mainland
#'   species without clades), \eqn{m_i = l_i + h_i},
#'   \eqn{o_i = \nu_i - \mathrm{own}_i} (without the branch being
#'   integrated) and \eqn{f(x) = \max(0, 1 - x/K_i)}: cladogenesis on
#'   \eqn{i} has weight \eqn{f(\nu_i + 1)} in the equations for E,
#'   \eqn{f(m_i + 1)} in those for F and \eqn{f(o_i + 1)} in those for DE
#'   and DM; colonisation of and dispersal into \eqn{i} have weight
#'   \eqn{f(\nu_i)} for E, \eqn{f(m_i)} for A and F and \eqn{f(o_i)} for DE
#'   and DM. A lineage thus counts itself when it speciates on its own island
#'   and not when it disperses to another island.
#'
#' @param dd Diversity-dependence data.
#' @param t Age.
#' @param l Contribution of the reconstructed lineages (per island).
#' @param own Contribution of the branch being integrated to \code{l}.
#' @return Matrix with one row per component and columns \code{E},
#'   \code{A}, \code{F} and \code{D}.
#' @keywords internal
DAISIE_A_dd_wmat <- function(dd, t, l, own = 0) {
  n <- length(dd$K)
  W <- matrix(1, 2 * n + 1, 4, dimnames = list(NULL, c("E", "A", "F", "D")))
  h <- DAISIE_A_dd_interp(dd$gt, dd$gh, t)
  p <- DAISIE_A_dd_interp(dd$gt, dd$gp, t)
  nu <- l + h + p
  m <- l + h
  o <- nu - own
  f <- function(x) ifelse(is.finite(dd$K), pmax(0, 1 - x / dd$K), 1)
  W[1 + seq_len(n), ] <- cbind(f(nu + 1), f(nu), f(m + 1), f(o + 1))
  W[1 + n + seq_len(n), ] <- cbind(f(nu), f(m), f(m), f(o))
  # weight of the recolonisation of the single island of the non-endemic
  # lineage (reset of its colonisation time; loss term of DM): f(o + 1)
  attr(W, "reset") <- f(o + 1)
  W
}

#' Linear interpolation on a grid of ages that may contain an age twice
#'
#' @param gt Increasing ages (an age may appear twice: the first copy holds
#'   the value on the younger side, the second the value on the older side).
#' @param M Matrix with one column per grid age.
#' @param t Age.
#' @return The interpolated column.
#' @keywords internal
DAISIE_A_dd_interp <- function(gt, M, t) {
  G <- length(gt)
  if (G == 1) return(M[, 1])
  g <- findInterval(t, gt, left.open = TRUE)
  g <- min(max(g, 1), G - 1)
  d <- gt[g + 1] - gt[g]
  a <- if (d > 0) (t - gt[g]) / d else 1
  a <- min(1, max(0, a))
  M[, g] + a * (M[, g + 1] - M[, g])
}

#' Contribution of the reconstructed lineages to nu at an age
#'
#' @param dd Diversity-dependence data of the likelihood context.
#' @param age Age; at a break (node or colonisation time) the value on the
#'   older side is returned.
#' @return Vector with one value per island.
#' @keywords internal
DAISIE_A_dd_l <- function(dd, age) {
  j <- findInterval(age, dd$lbreaks)
  dd$lvals[, j + 1]
}

#' Weights of the rate components at an age
#'
#' @param ctx Likelihood context.
#' @param age Age; at a node or colonisation time the diversity before the
#'   event (older side) is used.
#' @param block Block of equations: \code{"E"}, \code{"A"}, \code{"F"} or
#'   \code{"D"}.
#' @param own Contribution of the lineage concerned to the reconstructed
#'   lineages (per island).
#' @return \code{1} without diversity dependence, else one weight per rate
#'   component, see \code{\link{DAISIE_A_dd_wmat}()}.
#' @keywords internal
DAISIE_A_dd_weights <- function(ctx, age, block = "D", own = 0) {
  dd <- ctx$dd
  if (is.null(dd)) return(1)
  DAISIE_A_dd_wmat(dd, age, DAISIE_A_dd_l(dd, age), own = own)[, block]
}

#' Reconstructed lineages through time, by range of their sampled
#' descendants
#'
#' @param element Checked data of one mainland species.
#' @param range_index Named vector mapping range strings to indices.
#' @param island_age Age of the archipelago.
#' @return A list with \code{lbreaks} (increasing ages at which the
#'   reconstructed lineages change) and \code{Q}, a matrix with one column per
#'   interval between breaks (from the present), whose entry \eqn{R} is the
#'   sum over the reconstructed lineages alive in that interval of the
#'   fraction of their sampled descendants with range \eqn{R}.
#' @keywords internal
DAISIE_A_dd_tree <- function(element, range_index, island_age) {
  S <- length(range_index)
  segs <- list()
  clades <- element$clades
  for (cl in clades) {
    nd <- cl$nodes
    nr <- nrow(nd)
    counts <- matrix(0, nr, S)
    row_of <- match(nd$parent, nd$node)
    for (k in rev(seq_len(nr))) {
      if (!(nd$node[k] %in% nd$parent)) {
        counts[k, range_index[nd$range[k]]] <- 1
      }
      if (!is.na(row_of[k])) {
        counts[row_of[k], ] <- counts[row_of[k], ] + counts[k, ]
      }
    }
    for (k in seq_len(nr)) {
      to <- if (is.na(row_of[k])) cl$colonisation_time else nd$age[row_of[k]]
      segs[[length(segs) + 1]] <- list(from = nd$age[k], to = to,
                                       q = counts[k, ] / sum(counts[k, ]))
    }
  }
  if (length(clades) > 0 && !is.na(element$nonendemic_range)) {
    # non-endemic species with unknown colonisation time: it colonised after
    # the youngest clade; counted from the present to that colonisation time
    q <- rep(0, S)
    q[range_index[element$nonendemic_range]] <- 1
    t_young <- min(vapply(clades, function(cl) cl$colonisation_time,
                          numeric(1)))
    segs[[length(segs) + 1]] <- list(from = 0, to = t_young, q = q)
  }
  ends <- unlist(lapply(segs, function(s) c(s$from, s$to)))
  lbreaks <- sort(unique(ends[ends > 0 & ends < island_age]))
  bounds <- c(0, lbreaks, island_age)
  Q <- matrix(0, S, length(bounds) - 1)
  for (s in segs) {
    j <- which(bounds[-length(bounds)] >= s$from - 1e-12 &
                 bounds[-1] <= s$to + 1e-12)
    Q[, j] <- Q[, j] + s$q
  }
  list(lbreaks = lbreaks, Q = Q)
}

#' Mean-field contributions of the unobserved species
#'
#' @description Computes, for the current approximation \code{ctx$dd} of
#'   \eqn{\nu}, the expected numbers \eqn{h_R(t)} of species of the mainland
#'   species with range \eqn{R} at age \eqn{t} that do not belong to the
#'   reconstructed phylogenies and, for a mainland species without clades,
#'   the posterior probability \eqn{p_R(t)} that its non-endemic lineage is
#'   present with range \eqn{R}, on the grid \code{ctx$dd$gt}.
#'
#' @details The unobserved species follow the multi-type branching process
#'   conditioned on leaving no sampled descendants (rates of an event
#'   \eqn{R \to (a, b)} multiplied by \eqn{E_a E_b / E_R}, transitions
#'   \eqn{R \to R'} by \eqn{E_{R'}/E_R}, extinction divided by \eqn{E_R}).
#'   They are created by the reconstructed lineages (an event in which one
#'   daughter has no sampled descendants), by earlier colonisations without
#'   sampled descendants and, for a mainland species without clades, by its
#'   non-endemic lineage (posterior intensities from a forward-backward
#'   computation). The expected numbers are integrated as
#'   \eqn{G_R = h_R/E_R}, which removes the singular rates at the present.
#'
#' @param ctx Likelihood context with one power-series coefficient
#'   (\code{"rho"}-type initial conditions) and the current \code{dd}.
#' @param y0 Initial state (rows E, A, F) at the present.
#' @param tree Output of \code{\link{DAISIE_A_dd_tree}()}.
#' @param has_clades Logical: does the mainland species have clades?
#' @param t_col_old Colonisation time of the oldest clade (\code{Inf} if
#'   none).
#' @return List with matrices \code{H} and \code{p} (ranges x grid ages).
#' @keywords internal
DAISIE_A_dd_hidden <- function(ctx, y0, tree, has_clades, t_col_old) {
  dd <- ctx$dd
  gt <- dd$gt
  G <- length(gt)
  S <- length(ctx$ranges$sets)
  n <- ctx$pars$n_islands
  iE <- seq_len(S)
  iA <- S + 1
  iF <- S + 1 + seq_len(S)
  # pieces of the grid with a constant set of existing islands: a piece ends
  # at an island age (shared with the next piece) or at the first copy of a
  # submergence age (the next piece starts at the second copy)
  pieces <- list()
  start <- 1
  for (g in seq_len(G - 1) + 1) {
    if (gt[g] == gt[g - 1]) {
      pieces[[length(pieces) + 1]] <- start:(g - 1)
      start <- g
    } else if (gt[g] %in% ctx$breaks && g < G) {
      pieces[[length(pieces) + 1]] <- start:g
      start <- g
    }
  }
  pieces[[length(pieces) + 1]] <- start:G
  mid_l <- function(idx) {
    if (length(idx) < 2) return(matrix(0, n, 0))
    mids <- (gt[idx[-1]] + gt[idx[-length(idx)]]) / 2
    matrix(dd$lvals[, findInterval(mids, dd$lbreaks) + 1], nrow = n)
  }
  ptrs_at <- function(age) {
    comps <- DAISIE_A_dd_cpp_components(ctx, age)
    lapply(comps, function(st) st$cpp)
  }
  # backward pass on the grid (E, A, F); a repeated grid age is a submergence
  # age: first copy before, second copy after the jump (going back in time)
  B <- matrix(0, length(y0), G)
  y <- y0
  for (k in seq_along(pieces)) {
    idx <- pieces[[k]]
    if (k > 1 && min(idx) == max(pieces[[k - 1]]) + 1) {
      y <- as.vector(DAISIE_A_lik_submerge(matrix(y, ncol = 1), gt[min(idx)],
                                           ctx))
    }
    if (length(idx) == 1) {
      B[, idx] <- y
      next
    }
    res <- .Call("DAISIE_A_cpp_integrate_dd_grid", ptrs_at(mean(gt[idx])),
                 as.numeric(y), 1L, gt[idx], mid_l(idx),
                 "odeint::runge_kutta_cash_karp54", 1e-12, 1e-8,
                 as.numeric(dd$K), gt, dd$gh, dd$gp)
    B[, idx] <- res
    y <- res[, ncol(res)]
  }
  mid_Q <- function(idx) {
    if (length(idx) < 2) return(matrix(0, S, 0))
    mids <- (gt[idx[-1]] + gt[idx[-length(idx)]]) / 2
    matrix(tree$Q[, findInterval(mids, tree$lbreaks) + 1], nrow = S)
  }
  run <- function(mode, x0, src) {
    out <- matrix(0, length(x0), G)
    x <- x0
    for (k in rev(seq_along(pieces))) {
      idx <- pieces[[k]]
      if (k < length(pieces) && gt[max(idx)] == gt[min(pieces[[k + 1]])] &&
          max(idx) + 1 == min(pieces[[k + 1]])) {
        # submergence between this piece and the older one
        x <- DAISIE_A_dd_jump(x, gt[max(idx)], ctx, mode,
                              E_old = B[iE, max(idx) + 1],
                              E_new = B[iE, max(idx)])
      }
      if (length(idx) == 1) {
        out[, idx] <- x
        next
      }
      mids <- (gt[idx[-1]] + gt[idx[-length(idx)]]) / 2
      res <- .Call("DAISIE_A_cpp_dd_forward", ptrs_at(mean(gt[idx])),
                   as.integer(mode), gt[idx], B[iE, idx, drop = FALSE],
                   mid_l(idx), as.numeric(dd$K), gt, dd$gh, dd$gp,
                   mid_Q(idx),
                   as.integer(has_clades & mids > t_col_old),
                   if (is.null(src)) numeric(0) else src[, idx, drop = FALSE],
                   x, "odeint::runge_kutta_cash_karp54", 1e-10, 1e-7)
      out[, idx] <- res
      x <- res[, 1]
    }
    out
  }
  src_grid <- NULL
  pres_grid <- matrix(0, S, G)
  if (!has_clades) {
    # forward-backward for a mainland species without clades: posterior
    # presence of its non-endemic lineage and posterior intensity of the
    # creation of endemic species without sampled descendants
    U <- run(1L, c(1, rep(0, S)), NULL)
    src_grid <- matrix(0, S, G)
    for (g in seq_len(G)) {
      L <- sum(U[, g] * B[c(iA, iF), g])
      if (!is.finite(L) || L <= 0) next
      pres_grid[, g] <- U[1 + iE, g] * B[iF, g] / L
      t <- gt[g]
      W <- DAISIE_A_dd_wmat(dd, t, DAISIE_A_dd_l(dd, t))
      comps <- DAISIE_A_lik_get_components(ctx, t)
      E <- B[iE, g]
      A <- B[iA, g]
      for (c in seq_along(comps)) {
        st <- comps[[c]]
        wF <- W[c, "F"]
        if (wF <= 0) next
        if (length(st$bud_r) > 0) {
          v <- wF * st$bud_r * U[1 + st$bud_from, g] * B[iF[st$bud_p], g] / L
          src_grid[, g] <- src_grid[, g] + DAISIE_A_dd_tabulate(st$bud_x, v, S)
        }
        if (length(st$end_r) > 0) {
          v <- wF * st$end_r * U[1 + st$end_p, g] * A / L
          src_grid[, g] <- src_grid[, g] +
            DAISIE_A_dd_tabulate(st$end_a, v * E[st$end_b], S) +
            DAISIE_A_dd_tabulate(st$end_b, v * E[st$end_a], S)
        }
        if (st$laa > 0) {
          src_grid[, g] <- src_grid[, g] + st$laa * U[1 + iE, g] * A / L
        }
      }
    }
  }
  Gm <- run(0L, rep(0, S), src_grid)
  H <- B[iE, , drop = FALSE] * Gm
  H[H < 0] <- 0
  list(H = H, p = pres_grid)
}

#' Jump of the forward-in-time quantities of the mean-field approximation at
#' island submergence
#'
#' @param x State: \eqn{G = h/E} (\code{mode = 0}) or the adjoint
#'   \eqn{(U_A, U_F)} (\code{mode = 1}) just before the submergence.
#' @param age Submergence age.
#' @param ctx Likelihood context.
#' @param mode 0 or 1.
#' @param E_old,E_new E just before and just after the submergence (in
#'   forward time).
#' @return The state just after the submergence.
#' @keywords internal
DAISIE_A_dd_jump <- function(x, age, ctx, mode, E_old, E_new) {
  i <- which(abs(ctx$submergence_ages - age) <= 1e-12)
  if (mode == 0) {
    v <- x * E_old
  } else {
    absent <- x[1]
    v <- x[-1]
  }
  for (ii in i) {
    map <- ctx$sub_map[[ii]]
    v_old <- v
    v[which(!is.na(map))] <- 0
    for (s in which(!is.na(map))) {
      if (map[s] > 0) {
        v[map[s]] <- v[map[s]] + v_old[s]
      } else if (mode == 1) {
        absent <- absent + v_old[s]
      }
    }
  }
  if (mode == 0) {
    return(ifelse(E_new > 0, v / ifelse(E_new > 0, E_new, 1), 0))
  }
  c(absent, v)
}

#' Rate components with C++ structures
#'
#' @param ctx Likelihood context with diversity-dependence data.
#' @param age Age.
#' @return The rate components at \code{age}, each with an external pointer
#'   \code{cpp} to its C++ structures.
#' @keywords internal
DAISIE_A_dd_cpp_components <- function(ctx, age) {
  ctx$method <- "odeint::runge_kutta_cash_karp54"
  DAISIE_A_lik_get_components(ctx, age)
}

#' Sum weights by index
#'
#' @param idx Indices.
#' @param w Weights.
#' @param S Length of the result.
#' @return Vector of length \code{S}.
#' @keywords internal
DAISIE_A_dd_tabulate <- function(idx, w, S) {
  out <- rep(0, S)
  if (length(idx) == 0) return(out)
  s <- rowsum(w, idx)
  out[as.integer(rownames(s))] <- s[, 1]
  out
}

#' Mean-field diversity of one mainland species for finite carrying
#' capacities
#'
#' @param ctx Likelihood context (of the element).
#' @param element Checked data of the mainland species.
#' @param range_index Named vector mapping range strings to indices.
#' @param E0 Initial values of E at the present (one column,
#'   \code{"rho"}-type).
#' @param A0,F0 Initial values of A and F at the present (one column).
#' @param island_age Age of the archipelago.
#' @param control List with \code{grid} (number of equidistant grid ages),
#'   \code{tol} (convergence tolerance for \eqn{\nu}) and \code{maxit}.
#' @return The diversity-dependence data \code{dd} for the likelihood
#'   context: \code{comp_pars}, \code{K}, \code{lbreaks}, \code{lvals}
#'   (contribution of the reconstructed lineages per interval), \code{gt},
#'   \code{gh} and \code{gp} (grids of the contributions of the unobserved
#'   species and of the non-endemic lineage), and the number of iterations.
#' @keywords internal
DAISIE_A_dd_nu <- function(ctx, element, range_index, E0, A0, F0, island_age,
                           control) {
  pars <- ctx$pars
  n <- pars$n_islands
  occ <- matrix(t(vapply(ctx$ranges$sets, function(R) {
    as.numeric(seq_len(n) %in% R)
  }, numeric(n))), ncol = n)
  tree <- DAISIE_A_dd_tree(element, range_index, island_age)
  lvals <- t(occ) %*% tree$Q
  has_clades <- length(element$clades) > 0
  t_col_old <- if (has_clades) {
    max(vapply(element$clades, function(cl) cl$colonisation_time, numeric(1)))
  } else {
    Inf
  }
  sub_ages <- ctx$submergence_ages[is.finite(ctx$submergence_ages) &
                                     ctx$submergence_ages > 0]
  gt <- sort(c(unique(c(seq(0, island_age, length.out = control$grid),
                        tree$lbreaks, ctx$breaks)),
               sub_ages))
  zero <- matrix(0, n, length(gt))
  dd <- list(comp_pars = DAISIE_A_dd_components(pars),
             K = pars$K,
             lbreaks = tree$lbreaks,
             all_breaks = sort(unique(c(ctx$breaks, tree$lbreaks))),
             lvals = lvals,
             gt = gt,
             gh = zero,
             gp = zero)
  ctx1 <- ctx
  ctx1$poly <- DAISIE_A_poly()
  y0 <- rbind(E0, A0, F0)
  it <- 0
  repeat {
    it <- it + 1
    ctx1$dd <- dd
    res <- DAISIE_A_dd_hidden(ctx1, y0, tree, has_clades, t_col_old)
    gh_new <- t(occ) %*% res$H
    gp_new <- t(occ) %*% res$p
    change <- max(abs(c(gh_new - dd$gh, gp_new - dd$gp)))
    # damping after a few iterations, against oscillations
    damp <- if (it > 5) 0.5 else 1
    dd$gh <- dd$gh + damp * (gh_new - dd$gh)
    dd$gp <- dd$gp + damp * (gp_new - dd$gp)
    if (change < control$tol || it >= control$maxit) break
  }
  dd$iterations <- it
  dd$converged <- change < control$tol
  dd
}
