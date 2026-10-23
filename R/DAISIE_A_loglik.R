# Likelihood of the archipelago model DAISIE_A for infinite carrying
# capacities, generalising the D-E approach of DAISIE_DE.
#
# Quantities (all functions of the age t, integrated from the present t = 0
# backwards to the island age), one per allowed range R:
#   E_R   probability that an endemic lineage with range R at age t leaves no
#         species in the phylogeny at the present (it may leave missing
#         species, see below);
#   DE_R  probability (density) of the observed subtree below a branch, given
#         that the branch is an endemic lineage with range R at age t;
#   DM_R  the same, given that the branch is the non-endemic (mainland)
#         lineage with range R at age t, times the probability that the
#         future of the mainland species after it stops carrying the observed
#         subtree is consistent with the data;
#   F_R   probability that the future of the mainland species is consistent
#         with the data, given that it is present as a non-endemic lineage
#         with range R that does not carry observed species;
#   A     the same, given that the mainland species is absent from the
#         archipelago.
# With missing species and exact sampling, every quantity is a truncated
# multivariate power series in variables that count missing species per
# range (and of unknown range); it is stored as a vector of coefficients.

#' Rate structures for the archipelago likelihood
#'
#' @param pars A \code{DAISIE_A_pars} object.
#' @param exists 0/1 vector indicating which islands exist.
#' @param ranges Output of \code{\link{DAISIE_A_ranges}()}.
#' @return A list with transition matrices, out-rates and the lists of
#'   cladogenetic events.
#' @keywords internal
DAISIE_A_lik_structures <- function(pars, exists, ranges) {
  n <- pars$n_islands
  sets <- ranges$sets
  S <- length(sets)
  code_index <- integer(2^n - 1)
  code_index[ranges$codes] <- seq_len(S)
  idx <- function(set) code_index[sum(2^(set - 1))]
  single <- vapply(seq_len(n), function(i) idx(i), integer(1))
  D <- sweep(pars$D, 2, exists, "*")
  J <- sweep(pars$J, 2, exists, "*")
  col <- pars$gam * exists
  vic_split <- DAISIE_A_vic_per(pars) == "split"
  Qe <- matrix(0, S, S)
  Qm <- matrix(0, S, S)
  g <- rep(0, S)
  insitu <- rep(0, S)
  vic_tot <- rep(0, S)
  jump_tot <- rep(0, S)
  reset <- rep(0, S)
  tri <- matrix(numeric(0), ncol = 4)  # parent, daughter a, daughter b, rate
  bud <- matrix(numeric(0), ncol = 3)  # parent (stays), budded endemic, rate
  endt <- matrix(numeric(0), ncol = 4) # parent, endemic a, endemic b, rate
  for (s in seq_len(S)) {
    R <- sets[[s]]
    r <- length(R)
    for (j in setdiff(seq_len(n), R)) {
      Rj <- sort(c(R, j))
      if (DAISIE_A_range_allowed(Rj, pars)) {
        sj <- idx(Rj)
        dr <- sum(D[R, j])
        Qe[s, sj] <- Qe[s, sj] + dr
        Qm[s, sj] <- Qm[s, sj] + dr + col[j]
      }
      jr <- sum(J[R, j])
      if (jr > 0) {
        jump_tot[s] <- jump_tot[s] + jr
        tri <- rbind(tri, c(s, s, idx(j), jr))
        bud <- rbind(bud, c(s, idx(j), jr))
      }
    }
    if (r == 1) {
      g[s] <- pars$mu[R]
      reset[s] <- col[R]
    } else {
      for (i in R) {
        Ri <- setdiff(R, i)
        if (DAISIE_A_range_allowed(Ri, pars)) {
          si <- idx(Ri)
          Qe[s, si] <- Qe[s, si] + pars$mu[i]
          Qm[s, si] <- Qm[s, si] + pars$mu[i]
        }
      }
    }
    for (i in R) {
      insitu[s] <- insitu[s] + pars$lac[i]
      if (r == 1) {
        tri <- rbind(tri, c(s, s, s, pars$lac[i]))
        if (pars$clado_single == "split") {
          endt <- rbind(endt, c(s, s, s, pars$lac[i]))
        } else {
          bud <- rbind(bud, c(s, s, pars$lac[i]))
        }
      } else {
        tri <- rbind(tri, c(s, s, idx(i), pars$lac[i]))
        bud <- rbind(bud, c(s, idx(i), pars$lac[i]))
      }
    }
    if (r >= 2 && pars$lav > 0) {
      splits <- DAISIE_A_admissible_splits(R, pars)
      ns <- length(splits)
      vt <- if (vic_split) pars$lav * ns else pars$lav
      vic_tot[s] <- vt
      for (part_a in splits) {
        a <- idx(sort(part_a))
        b <- idx(sort(setdiff(R, part_a)))
        tri <- rbind(tri, c(s, a, b, vt / ns))
        endt <- rbind(endt, c(s, a, b, vt / ns))
      }
    }
  }
  # drop events with zero rate
  tri <- tri[tri[, 4] > 0, , drop = FALSE]
  bud <- bud[bud[, 3] > 0, , drop = FALSE]
  endt <- endt[endt[, 4] > 0, , drop = FALSE]
  agg <- function(p) {
    m <- matrix(0, S, length(p))
    m[cbind(p, seq_along(p))] <- 1
    m
  }
  omega_E <- rowSums(Qe) + g + insitu + vic_tot + jump_tot
  omega_F <- rowSums(Qm) + g + insitu + vic_tot + jump_tot + pars$laa
  list(
    S = S,
    single = single,
    col = col,
    Qe = Qe,
    Qm = Qm,
    g = g,
    omega_E = omega_E,
    omega_F = omega_F,
    omega_M = omega_F + reset,
    laa = pars$laa,
    tri_p = tri[, 1], tri_a = tri[, 2], tri_b = tri[, 3], tri_r = tri[, 4],
    tri_same = tri[, 2] == tri[, 3],
    A_tri = agg(tri[, 1]),
    bud_p = bud[, 1], bud_x = bud[, 2], bud_r = bud[, 3],
    A_bud = agg(bud[, 1]),
    end_p = endt[, 1], end_a = endt[, 2], end_b = endt[, 3],
    end_r = endt[, 4],
    end_same = endt[, 2] == endt[, 3],
    A_end = agg(endt[, 1])
  )
}

#' Truncated multivariate power series used for missing species
#'
#' @param degrees Integer vector with the maximum degree of every variable
#'   (\code{integer(0)} for no variables).
#' @return A list with the number of coefficients \code{P}, the multi-indices
#'   \code{K}, a function \code{index()} giving the position of a
#'   multi-index, and a function \code{mul(U, V)} multiplying the series in
#'   the rows of \code{U} and \code{V} (matrices with \code{P} columns;
#'   \code{V} may have a single row).
#' @keywords internal
DAISIE_A_poly <- function(degrees = integer(0)) {
  degrees <- as.integer(degrees)
  if (length(degrees) == 0) {
    return(list(P = 1L, K = matrix(0L, 1, 0), degrees = degrees,
                index = function(k) 1L,
                mul = function(U, V) {
                  if (nrow(V) == 1 && nrow(U) > 1) V <- V[rep(1, nrow(U)), , drop = FALSE]
                  U * V
                }))
  }
  K <- as.matrix(expand.grid(lapply(degrees, function(d) 0:d)))
  dimnames(K) <- NULL
  P <- nrow(K)
  radix <- cumprod(c(1, degrees + 1))[seq_along(degrees)]
  index <- function(k) as.integer(1 + sum(k * radix))
  pairs <- which(matrix(TRUE, P, P), arr.ind = TRUE)
  ksum <- K[pairs[, 1], , drop = FALSE] + K[pairs[, 2], , drop = FALSE]
  ok <- apply(ksum, 1, function(k) all(k <= degrees))
  pairs <- pairs[ok, , drop = FALSE]
  ksum <- ksum[ok, , drop = FALSE]
  ic <- as.integer(1 + ksum %*% radix)
  Cmat <- matrix(0, nrow(pairs), P)
  Cmat[cbind(seq_len(nrow(pairs)), ic)] <- 1
  ia <- pairs[, 1]
  ib <- pairs[, 2]
  mul <- function(U, V) {
    if (nrow(V) == 1 && nrow(U) > 1) V <- V[rep(1, nrow(U)), , drop = FALSE]
    (U[, ia, drop = FALSE] * V[, ib, drop = FALSE]) %*% Cmat
  }
  list(P = P, K = K, degrees = degrees, index = index, mul = mul)
}

#' Right-hand side of the archipelago likelihood equations
#'
#' @param t Age (not used: the rates are constant within a segment).
#' @param y State vector: the matrix with rows E (S), A (1), F (S) and, if
#'   present, DE (S) and DM (S), and one column per power-series coefficient.
#' @param parms List with the rate structures \code{st} and the power-series
#'   algebra \code{poly}.
#' @return List with the derivative, as required by \code{deSolve::ode()}.
#' @keywords internal
DAISIE_A_lik_rhs <- function(t, y, parms) {
  st <- parms$st
  mul <- parms$poly$mul
  P <- parms$poly$P
  S <- st$S
  m <- matrix(y, ncol = P)
  has_D <- nrow(m) > 2 * S + 1
  E <- m[seq_len(S), , drop = FALSE]
  A <- m[S + 1, , drop = FALSE]
  Fm <- m[S + 1 + seq_len(S), , drop = FALSE]
  one <- matrix(c(1, rep(0, P - 1)), nrow = 1)
  # endemic lineages without observed descendants
  Eab <- mul(E[st$tri_a, , drop = FALSE], E[st$tri_b, , drop = FALSE])
  dE <- -st$omega_E * E + st$Qe %*% E + st$g %*% one +
    st$A_tri %*% (st$tri_r * Eab)
  # mainland species absent
  dA <- -sum(st$col) * A + st$col %*% Fm[st$single, , drop = FALSE]
  # mainland species present as a non-endemic not carrying observed species
  FE <- mul(Fm[st$bud_p, , drop = FALSE], E[st$bud_x, , drop = FALSE])
  endE <- st$laa * E +
    st$A_end %*% (st$end_r * mul(E[st$end_a, , drop = FALSE],
                                 E[st$end_b, , drop = FALSE]))
  dF <- -st$omega_F * Fm + st$Qm %*% Fm + st$g %*% A +
    st$A_bud %*% (st$bud_r * FE) + mul(endE, A)
  if (!has_D) {
    return(list(as.vector(rbind(dE, dA, dF))))
  }
  DE <- m[2 * S + 1 + seq_len(S), , drop = FALSE]
  DM <- m[3 * S + 1 + seq_len(S), , drop = FALSE]
  DEab <- mul(DE[st$tri_a, , drop = FALSE], E[st$tri_b, , drop = FALSE]) +
    mul(E[st$tri_a, , drop = FALSE], DE[st$tri_b, , drop = FALSE])
  dDE <- -st$omega_E * DE + st$Qe %*% DE + st$A_tri %*% (st$tri_r * DEab)
  DMbud <- mul(DM[st$bud_p, , drop = FALSE], E[st$bud_x, , drop = FALSE]) +
    mul(Fm[st$bud_p, , drop = FALSE], DE[st$bud_x, , drop = FALSE])
  DEend <- mul(DE[st$end_a, , drop = FALSE], E[st$end_b, , drop = FALSE]) +
    mul(E[st$end_a, , drop = FALSE], DE[st$end_b, , drop = FALSE])
  endD <- st$laa * DE + st$A_end %*% (st$end_r * DEend)
  dDM <- -st$omega_M * DM + st$Qm %*% DM + st$A_bud %*% (st$bud_r * DMbud) +
    mul(endD, A)
  list(as.vector(rbind(dE, dA, dF, dDE, dDM)))
}

#' Integrate the archipelago likelihood equations between two ages
#'
#' @param y State matrix (rows E, A, F[, DE, DM]; one column per
#'   coefficient).
#' @param from,to Ages (\code{from < to}).
#' @param ctx Likelihood context (see \code{\link{DAISIE_A_lik_context}()}).
#' @return The state matrix at age \code{to}.
#' @keywords internal
DAISIE_A_lik_integrate <- function(y, from, to, ctx) {
  if (to <= from) return(y)
  breaks <- ctx$breaks[ctx$breaks > from & ctx$breaks < to]
  bounds <- c(from, breaks, to)
  nr <- nrow(y)
  for (k in seq_len(length(bounds) - 1)) {
    st <- DAISIE_A_lik_get_structures(ctx, (bounds[k] + bounds[k + 1]) / 2)
    sol <- deSolve::ode(
      y = as.vector(y),
      times = bounds[k:(k + 1)],
      func = DAISIE_A_lik_rhs,
      parms = list(st = st, poly = ctx$poly),
      method = ctx$method,
      rtol = ctx$rtol,
      atol = ctx$atol
    )
    if (nrow(sol) < 2 || any(!is.finite(sol[2, -1]))) {
      stop("Integration of the likelihood equations failed.")
    }
    y <- matrix(sol[2, -1], nrow = nr)
  }
  return(y)
}

#' Rate structures at a given age (cached per set of existing islands)
#'
#' @param ctx Likelihood context.
#' @param age Age (time before present).
#' @return Rate structures, see \code{\link{DAISIE_A_lik_structures}()}.
#' @keywords internal
DAISIE_A_lik_get_structures <- function(ctx, age) {
  exists <- as.numeric(age <= ctx$island_ages + 1e-12)
  key <- paste(exists, collapse = "")
  st <- ctx$cache[[key]]
  if (is.null(st)) {
    st <- DAISIE_A_lik_structures(ctx$pars, exists, ctx$ranges)
    assign(key, st, envir = ctx$cache)
  }
  return(st)
}

#' Context of a likelihood computation
#'
#' @param pars A \code{DAISIE_A_pars} object.
#' @param island_age Age of the archipelago (oldest island).
#' @param poly Power-series algebra, see \code{\link{DAISIE_A_poly}()}.
#' @param methode,rtol,atol ODE solver settings.
#' @param cache Environment for caching rate structures.
#' @return A list.
#' @keywords internal
DAISIE_A_lik_context <- function(pars, island_age, poly, methode, rtol, atol,
                                 cache = new.env()) {
  island_ages <- if (is.null(pars$island_ages)) {
    rep(island_age, pars$n_islands)
  } else {
    pars$island_ages
  }
  if (any(island_ages > island_age + 1e-12)) {
    stop("Island ages cannot be larger than the age of the archipelago.")
  }
  list(
    pars = pars,
    ranges = DAISIE_A_ranges(pars),
    island_ages = island_ages,
    breaks = sort(unique(island_ages[island_ages < island_age])),
    poly = poly,
    method = methode,
    rtol = rtol,
    atol = atol,
    cache = cache
  )
}

#' Log-likelihood of the data of one mainland species under DAISIE_A
#'
#' @description Computes the log-likelihood of the colonisation times,
#'   reconstructed phylogenies, tip ranges and endemicity of all descendants
#'   on the archipelago of one mainland species, for infinite carrying
#'   capacities. See \code{\link{DAISIE_A_loglik}()} for the data format.
#'
#' @param pars A \code{DAISIE_A_pars} object (K is ignored).
#' @param element Data of one mainland species (an element of an archipelago
#'   datalist), or \code{NULL} for a mainland species without descendants on
#'   the archipelago.
#' @param island_age Age of the archipelago.
#' @param sampling \code{"exact"} (default): exact likelihood of the observed
#'   numbers of missing species per range; \code{"rho"}: every species with
#'   range \eqn{R} is included in the phylogeny independently with
#'   probability \eqn{\rho_R}, estimated as the observed fraction.
#' @param methode ODE solver passed to \code{deSolve::ode()}.
#' @param rtol,atol Relative and absolute tolerance of the ODE solver.
#' @param ctx_cache Optional environment used to cache rate structures
#'   between calls with the same parameters.
#'
#' @return The log-likelihood.
#' @family DAISIE_A
#' @export
DAISIE_A_loglik_element <- function(pars,
                                    element,
                                    island_age,
                                    sampling = "exact",
                                    methode = "lsoda",
                                    rtol = 1e-10,
                                    atol = 1e-14,
                                    ctx_cache = new.env()) {
  sampling <- match.arg(sampling, c("exact", "rho"))
  n <- pars$n_islands
  ranges <- DAISIE_A_ranges(pars)
  S <- length(ranges$sets)
  range_index <- stats::setNames(seq_len(S), ranges$strings)
  if (is.null(element)) {
    element <- list(clades = list(), nonendemic_range = NA, missing = c())
  }
  element <- DAISIE_A_check_element(element, n)
  clades <- element$clades
  # sampled endemic tips per range
  tip_ranges <- unlist(lapply(clades, function(cl) {
    tips <- cl$nodes[!is.na(cl$nodes$range) & cl$nodes$type != "I", ]
    tips$range
  }))
  if (any(!(c(tip_ranges, element$missing[!is.na(element$missing)]) %in%
            ranges$strings))) {
    stop("A range in the data is not an allowed range under 'pars'.")
  }
  S_R <- tabulate(range_index[tip_ranges], nbins = S)
  miss <- element$missing
  m_R <- tabulate(range_index[miss[!is.na(miss)]], nbins = S)
  m_u <- sum(is.na(miss))
  # initial values of E and the polynomial algebra
  if (sampling == "exact") {
    miss_ranges <- which(m_R > 0)
    degrees <- c(m_R[miss_ranges], if (m_u > 0) m_u)
    poly <- DAISIE_A_poly(degrees)
    E0 <- matrix(0, S, poly$P)
    nv <- length(degrees)
    unit <- function(v) {
      k <- rep(0L, nv)
      k[v] <- 1L
      poly$index(k)
    }
    for (v in seq_along(miss_ranges)) {
      E0[miss_ranges[v], unit(v)] <- 1
    }
    if (m_u > 0) E0[, unit(nv)] <- E0[, unit(nv)] + 1
    target <- poly$index(c(m_R[miss_ranges], if (m_u > 0) m_u))
    log_tip <- 0
  } else {
    poly <- DAISIE_A_poly()
    m_alloc <- m_R
    if (m_u > 0) {
      w <- if (sum(S_R) > 0) S_R / sum(S_R) else rep(1 / S, S)
      m_alloc <- m_alloc + m_u * w
    }
    rho <- ifelse(S_R + m_alloc > 0, S_R / (S_R + m_alloc), 1)
    E0 <- matrix(1 - rho, S, 1)
    target <- 1L
    log_tip <- sum(S_R[S_R > 0] * log(rho[S_R > 0]))
  }
  P <- poly$P
  ctx <- DAISIE_A_lik_context(pars, island_age, poly, methode, rtol, atol,
                              cache = ctx_cache)
  zero <- matrix(0, 1, P)
  zeroS <- matrix(0, S, P)
  one <- matrix(c(1, rep(0, P - 1)), 1, P)
  logscale <- 0
  # initial conditions of the mainland species at the present, for the
  # youngest clade (post-phase)
  has_I <- any(vapply(clades, function(cl) any(cl$nodes$type %in% "I"),
                      logical(1)))
  A0 <- one
  F0 <- zeroS
  if (has_I) {
    A0 <- zero
  } else if (!is.na(element$nonendemic_range)) {
    A0 <- zero
    F0[range_index[element$nonendemic_range], 1] <- 1
  }
  pre <- NULL
  t_next <- Inf
  E_at <- NULL
  nc <- length(clades)
  if (nc > 0) {
    for (c in rev(seq_len(nc))) {
      env0 <- if (c == nc) rbind(E0, A0, F0) else rbind(E0, zero, zeroS)
      res <- DAISIE_A_clade_lik(
        clade = clades[[c]],
        env0 = env0,
        t_next = t_next,
        pre_next = pre,
        ctx = ctx,
        range_index = range_index
      )
      if (!is.finite(res$logscale)) return(-Inf)
      pre <- res$pre
      t_next <- clades[[c]]$colonisation_time
      logscale <- logscale + res$logscale
    }
    y <- pre
    from <- t_next
  } else {
    y <- rbind(E0, A0, F0)
    from <- 0
  }
  y <- DAISIE_A_lik_integrate(y, from, island_age, ctx)
  lik <- y[S + 1, target]
  if (!is.finite(lik) || lik <= 0) return(-Inf)
  loglik <- log(lik) + logscale + log_tip
  if (sampling == "exact") {
    S_tot <- sum(S_R)
    m_tot <- sum(m_R) + m_u
    loglik <- loglik - lchoose(S_tot + m_tot, S_tot)
  }
  return(loglik)
}

#' Likelihood contribution of one clade (colonisation event)
#'
#' @param clade A clade: list with \code{colonisation_time} and \code{nodes}.
#' @param env0 Environment (rows E, A, F) at the present for this clade.
#' @param t_next Colonisation time of the next younger clade (\code{Inf} if
#'   none).
#' @param pre_next Pre-phase state (rows E, A, F) of the next younger clade at
#'   \code{t_next}.
#' @param ctx Likelihood context.
#' @param range_index Named vector mapping range strings to indices.
#' @return List with \code{pre}, the state (rows E, A, F) at the colonisation
#'   time of this clade (the pre-phase initial condition), and the log of the
#'   scaling factors (\code{-Inf} if the data have probability zero).
#' @keywords internal
DAISIE_A_clade_lik <- function(clade, env0, t_next, pre_next, ctx,
                               range_index) {
  S <- length(range_index)
  P <- ctx$poly$P
  mul <- ctx$poly$mul
  nodes <- clade$nodes
  t_col <- clade$colonisation_time
  logscale <- 0
  zeroS <- matrix(0, S, P)
  iE <- seq_len(S)
  iA <- S + 1
  iF <- S + 1 + seq_len(S)
  iDE <- 2 * S + 1 + seq_len(S)
  iDM <- 3 * S + 1 + seq_len(S)
  # integrate a branch, injecting the pre-phase of the younger clade
  branch <- function(y, from, to) {
    if (t_next > from && t_next < to) {
      y <- DAISIE_A_lik_integrate(y, from, t_next, ctx)
      y[c(iA, iF), ] <- pre_next[c(iA, iF), ]
      y[iE, ] <- pre_next[iE, ]
      y <- DAISIE_A_lik_integrate(y, t_next, to, ctx)
    } else {
      y <- DAISIE_A_lik_integrate(y, from, to, ctx)
    }
    y
  }
  children <- split(nodes$node[nodes$parent != 0],
                    nodes$parent[nodes$parent != 0])
  state <- list()
  # post-order: nodes are numbered so that parents precede children
  for (k in rev(seq_len(nrow(nodes)))) {
    id <- nodes$node[k]
    age <- nodes$age[k]
    ch <- children[[as.character(id)]]
    if (is.null(ch)) {
      # tip
      y <- rbind(env0, zeroS, zeroS)
      ri <- range_index[nodes$range[k]]
      if (nodes$type[k] == "I") {
        y[iDM[ri], 1] <- 1
      } else {
        y[iDE[ri], 1] <- 1
      }
      if (age > 0) y <- branch(y, 0, age)
    } else {
      if (length(ch) != 2) stop("The phylogeny must be bifurcating.")
      yl <- branch(state[[as.character(ch[1])]],
                   nodes$age[nodes$node == ch[1]], age)
      yr <- branch(state[[as.character(ch[2])]],
                   nodes$age[nodes$node == ch[2]], age)
      st <- DAISIE_A_lik_get_structures(ctx, age)
      y <- yl
      y[iE, ] <- (yl[iE, ] + yr[iE, ]) / 2
      DEl <- yl[iDE, , drop = FALSE]
      DEr <- yr[iDE, , drop = FALSE]
      DMl <- yl[iDM, , drop = FALSE]
      DMr <- yr[iDM, , drop = FALSE]
      pair <- function(a, b, same) {
        w <- mul(DEl[a, , drop = FALSE], DEr[b, , drop = FALSE])
        w2 <- mul(DEl[b, , drop = FALSE], DEr[a, , drop = FALSE])
        ifelse(same, 1, 0.5) * w + ifelse(same, 0, 0.5) * w2
      }
      DEn <- st$A_tri %*% (st$tri_r * pair(st$tri_a, st$tri_b, st$tri_same))
      DMb <- 0.5 * (mul(DMl[st$bud_p, , drop = FALSE],
                        DEr[st$bud_x, , drop = FALSE]) +
                      mul(DEl[st$bud_x, , drop = FALSE],
                          DMr[st$bud_p, , drop = FALSE]))
      DMend <- st$A_end %*% (st$end_r * pair(st$end_a, st$end_b, st$end_same))
      DMn <- st$A_bud %*% (st$bud_r * DMb) +
        mul(DMend, y[iA, , drop = FALSE])
      sc <- max(abs(c(DEn, DMn)))
      if (!is.finite(sc) || sc <= 0) {
        return(list(pre = NULL, logscale = -Inf))
      }
      y[iDE, ] <- DEn / sc
      y[iDM, ] <- DMn / sc
      logscale <- logscale + log(sc)
      state[[as.character(ch[1])]] <- NULL
      state[[as.character(ch[2])]] <- NULL
    }
    state[[as.character(id)]] <- y
  }
  crown <- nodes$node[nodes$parent == 0]
  y <- branch(state[[as.character(crown)]], nodes$age[nodes$parent == 0],
              t_col)
  # colonisation at t_col: the stem is a non-endemic lineage on one island
  st <- DAISIE_A_lik_get_structures(ctx, t_col)
  DMs <- y[iDM[st$single], , drop = FALSE]
  pre <- y[c(iE, iA, iF), , drop = FALSE]
  pre[iA, ] <- st$col %*% DMs
  pre[iF, ] <- 0
  pre[iF[st$single], ] <- st$col * DMs
  sc <- max(abs(pre[c(iA, iF), ]))
  if (!is.finite(sc) || sc <= 0) {
    return(list(pre = NULL, logscale = -Inf))
  }
  pre[c(iA, iF), ] <- pre[c(iA, iF), ] / sc
  logscale <- logscale + log(sc)
  list(pre = pre, logscale = logscale)
}

#' Check and standardise the data of one mainland species
#'
#' @param element Data element.
#' @param n_islands Number of islands.
#' @return The standardised element (clades sorted from old to young).
#' @keywords internal
DAISIE_A_check_element <- function(element, n_islands) {
  if (is.null(element$clades)) element$clades <- list()
  if (is.null(element$nonendemic_range)) element$nonendemic_range <- NA
  if (is.null(element$missing)) element$missing <- character(0)
  clades <- element$clades
  if (length(clades) > 0) {
    col_times <- vapply(clades, function(cl) cl$colonisation_time,
                        numeric(1))
    clades <- clades[order(col_times, decreasing = TRUE)]
    for (k in seq_along(clades)) {
      nd <- clades[[k]]$nodes
      if (is.null(nd) || !all(c("node", "parent", "age", "range", "type") %in%
                              names(nd))) {
        stop("Every clade needs a node table with columns node, parent, age, ",
             "range and type.")
      }
      nd <- nd[order(nd$node), ]
      if (any(nd$parent >= nd$node)) {
        stop("Nodes must be numbered such that parents precede children.")
      }
      if (sum(nd$parent == 0) != 1) stop("A clade needs exactly one root.")
      if (clades[[k]]$colonisation_time < max(nd$age)) {
        stop("The colonisation time must be older than the crown.")
      }
      tips <- !(nd$node %in% nd$parent)
      if (any(is.na(nd$range[tips])) ||
          any(nchar(nd$range[tips]) != n_islands)) {
        stop("Every tip needs a range string of length n_islands.")
      }
      if (any(!(nd$type[tips] %in% c("I", "E")))) {
        stop("Tip types must be 'I' (non-endemic) or 'E' (endemic).")
      }
      nd$range[!tips] <- NA
      nd$type[!tips] <- NA
      if (sum(nd$type %in% "I") > 1) {
        stop("A clade can contain at most one non-endemic species.")
      }
      if (any(nd$type %in% "I") && k < length(clades)) {
        stop("Only the youngest clade can contain the non-endemic species.")
      }
      clades[[k]]$nodes <- nd
    }
    if (any(vapply(clades, function(cl) any(cl$nodes$type %in% "I"),
                   logical(1))) && !is.na(element$nonendemic_range)) {
      stop("The non-endemic species is either a tip of the youngest clade ",
           "or present with unknown colonisation time, not both.")
    }
  }
  element$clades <- clades
  return(element)
}

#' Log-likelihood of the archipelago model for infinite carrying capacity
#'
#' @description
#' Computes the log-likelihood of archipelago data under DAISIE_A with
#' \eqn{K = \infty}, using the D-E approach of DAISIE_DE generalised to ranges
#' (see the report for the equations). The data of every mainland species
#' consist of the colonisation times, the reconstructed phylogenies with tip
#' ranges and endemicity, and the ranges of missing species (species known to
#' be present but absent from the phylogeny). For a single island the
#' log-likelihood equals that of \code{\link{DAISIE_DE_loglik_CS}()}.
#'
#' Recolonisation convention: the colonisation time of a non-endemic lineage
#' is only reset when the mainland species recolonises the single island on
#' which it occurs; colonisation of another island expands its range.
#'
#' @param pars A \code{DAISIE_A_pars} object (K is ignored).
#' @param datalist Archipelago datalist: element \code{[[1]]} contains
#'   \code{island_age} and \code{not_present} (number of mainland species
#'   without descendants on the archipelago); elements \code{[[2]]},
#'   \code{[[3]]}, ... describe one mainland species each, with
#'   \itemize{
#'     \item \code{clades}: list of colonisation events with surviving
#'       descendants, each a list with \code{colonisation_time} (age) and
#'       \code{nodes}, a data frame with one row per node of the
#'       reconstructed tree: \code{node} (integer ID, parents before
#'       children), \code{parent} (0 for the root), \code{age},
#'       \code{range} (range string, e.g. \code{"0110"}, for tips) and
#'       \code{type} (\code{"E"} endemic or \code{"I"} non-endemic, for tips).
#'       At most one tip, in the youngest clade, can be non-endemic;
#'     \item \code{nonendemic_range} (optional): range of the non-endemic
#'       species if it is present but its colonisation time is unknown (and
#'       it is not part of a clade); \code{NA} otherwise;
#'     \item \code{missing} (optional): character vector with the ranges of
#'       endemic species of this mainland species that are missing from the
#'       phylogenies (\code{NA} for an unknown range).
#'   }
#'   See \code{\link{DAISIE_A_sim_to_datalist}()}.
#' @param sampling \code{"exact"} or \code{"rho"}, see
#'   \code{\link{DAISIE_A_loglik_element}()}.
#' @param cond Conditioning: 0 none, 1 on at least one colonisation with
#'   surviving descendants (as in DAISIE).
#' @param methode,rtol,atol ODE solver settings.
#' @param verbose Logical; print the log-likelihood per mainland species.
#'
#' @return The log-likelihood.
#' @family DAISIE_A
#' @export
#' @examples
#' pars <- DAISIE_A_create_pars(n_islands = 2, lac = 0.5, mu = 0.3,
#'   gam = 0.05, laa = 0.5, lav = 0.3, disp = 0.2)
#' set.seed(1)
#' sim <- DAISIE_A_sim(time = 2, M = 20, pars = pars, replicates = 1,
#'   verbose = FALSE)
#' datalist <- DAISIE_A_sim_to_datalist(sim[[1]])
#' DAISIE_A_loglik(pars, datalist)
DAISIE_A_loglik <- function(pars,
                            datalist,
                            sampling = "exact",
                            cond = 0,
                            methode = "lsoda",
                            rtol = 1e-10,
                            atol = 1e-14,
                            verbose = FALSE) {
  if (!DAISIE_A_are_pars(pars)) {
    stop("'pars' must be created with DAISIE_A_create_pars().")
  }
  island_age <- datalist[[1]]$island_age
  not_present <- datalist[[1]]$not_present
  cache <- new.env()
  logp0 <- DAISIE_A_loglik_element(pars, NULL, island_age,
                                   sampling = sampling, methode = methode,
                                   rtol = rtol, atol = atol,
                                   ctx_cache = cache)
  loglik <- not_present * logp0
  numimm <- not_present + length(datalist) - 1
  if (cond == 1) {
    loglik <- loglik - log(1 - exp(numimm * logp0))
  }
  for (k in seq_len(length(datalist) - 1) + 1) {
    lk <- DAISIE_A_loglik_element(pars, datalist[[k]], island_age,
                                  sampling = sampling, methode = methode,
                                  rtol = rtol, atol = atol,
                                  ctx_cache = cache)
    if (verbose) message("Mainland species ", k - 1, ": loglik = ", lk)
    loglik <- loglik + lk
  }
  return(loglik)
}
