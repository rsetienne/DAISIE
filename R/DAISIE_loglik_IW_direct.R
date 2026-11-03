# The IW likelihood with the right-hand side and the event step evaluated on the fly (rhs = "direct" of DAISIE_loglik_IW()).
#
# DAISIE_loglik_IW() with rhs = "tables" builds, for every ODE segment, eight coefficient arrays, six helper arrays and the
# matrix ki (all of the size of the state space or its square, see DAISIE_IW_pars(), nndivdep() and create_l0ki()), and the event
# between two segments works on further arrays of the size of the state space.  Here nothing of that size is created in R:
# the right-hand side and the event step are in C++ (src/DAISIE_IW_direct.cpp) and compute the coefficients from the state
# itself.  R only holds the vector of probabilities, and passes a few scalars.  The state (n1, n2, configuration a of the open
# clades) is stored at 1 + n1 + lxm * (n2 + lxe * a), the same numbering as the array with dim c(lxm, lxe, sysdim) of the
# tables version.
#
# rhs = "compact" is rhs = "direct" with a smaller state space.  The states with i + j > ceiling(K') - k (k = number of events
# so far) are exactly 0 for a finite K', because every transition that increases i + j has the factor max(0, 1 - (i + j + k) / K')
# (see the comments of src/DAISIE_IW_direct.cpp), so they are not stored.  Per configuration the states are then stored per j
# with i = 0..len[j]-1, at 1 + a * P + off[j] + i (DAISIE_IW_layout()); with K' = Inf, or if K' is larger than the box, this is the
# numbering of "direct" and of the tables.

# The largest number of states that rhs = "direct" supports: R's integer limit, 2^31 - 1.
DAISIE_IW_max_states <- function() {
  .Machine$integer.max
}

# Error with a clear message if the state space (lxm x lxe x sysdim states, or P x sysdim if P states per configuration are stored,
# as in rhs = "compact") has more than 2^31 - 1 states.  The product is computed in double precision, so it cannot overflow.
DAISIE_IW_check_dim <- function(lxm, lxe, sysdim, P = as.numeric(lxm) * as.numeric(lxe)) {
  nstates <- as.numeric(P) * as.numeric(sysdim)
  if (!is.finite(nstates) || nstates > DAISIE_IW_max_states()) {
    what <- if (P == as.numeric(lxm) * as.numeric(lxe)) {
      paste0("lxm x lxe x 2^s = ", format(lxm, scientific = FALSE), " x ", format(lxe, scientific = FALSE), " x ",
             format(sysdim, scientific = FALSE))
    } else {
      paste0("(states per configuration) x 2^s = ", format(P, scientific = FALSE), " x ", format(sysdim, scientific = FALSE))
    }
    stop("The state space of the island-wide likelihood has ", what, " = ", format(nstates, scientific = TRUE, digits = 4),
         " states, more than the 2^31 - 1 = ", DAISIE_IW_max_states(), " states that rhs = \"direct\" and rhs = \"compact\" can index. ",
         "Reduce the size of the state space (pars2[1], or the number of lineages that are expanded at the same time).",
         call. = FALSE)
  }
  invisible(nstates)
}

# The layout of the states of one configuration of the open clades when only the states with i + j <= tri are stored (tri = Inf:
# all states of the box): row j holds i = 0..len[j]-1, starting at off[j] (0-based); P states per configuration.  The state
# (i, j, a) is at 1 + a * P + off[j] + i.  The same as the C++ layout (daisie_iw_layout()).
DAISIE_IW_layout <- function(lxm, lxe, tri = Inf) {
  len <- pmin(lxm, tri - (0:(lxe - 1)) + 1)
  len <- len[len > 0]                           # len does not increase with j, so these are the first rows
  nrow <- length(len)
  list(len = len, off = if (nrow > 0) c(0, cumsum(len))[seq_len(nrow)] else numeric(0), P = sum(len))
}

# The bound tri of the segment with k events so far: the states with i + j > ceiling(K') - k have probability 0 (K' finite).
DAISIE_IW_tri <- function(Kprime, k, compact) {
  if (compact && is.finite(Kprime)) ceiling(Kprime) - k else Inf
}

# The number of open clades at the start of each segment k = 0..nbrts-2 (as DAISIE_IW_max_open_clades(), for every segment).
DAISIE_IW_segments <- function(event, nbrts) {
  s <- numeric(nbrts - 1)
  open <- 0
  for (k in 0:(nbrts - 2)) {
    s[k + 1] <- open
    if (k < nbrts - 2) {
      if (event[k + 2] == 1) {
        open <- open + 1
      } else if (event[k + 2] == 2) {
        open <- open - 1
      }
    }
  }
  s
}

# The largest number of clades that are open (colonised but without a first speciation) at the same time, s; the state space
# has at most lxm x lxe x 2^s states.  Follows the events of the loop of DAISIE_loglik_IW(): a colonisation (event 1) opens a
# clade, the first speciation of a clade (event 2) closes it.
DAISIE_IW_max_open_clades <- function(event, nbrts) {
  s <- 0
  smax <- 0
  if (nbrts >= 3) {
    for (k in 0:(nbrts - 3)) {
      if (event[k + 2] == 1) {
        s <- s + 1
      } else if (event[k + 2] == 2) {
        s <- s - 1
      }
      smax <- max(smax, s)
    }
  }
  return(smax)
}

# As checkprobs2(), with the work done in C++: values below 0 are set to 0; the vector is normalised and the log of its sum is
# added to the loglik, or the loglik becomes -Inf if the vector has NaN's or a sum that is not positive.
DAISIE_IW_check_probs <- function(probs, loglik, verbose) {
  res <- .Call("daisie_iw_normalize", probs)
  if (res$status == 0) {
    loglik <- loglik + log(res$sum)
  } else {
    loglik <- -Inf
    if (verbose) {
      message("Numerical issues encountered.")
    }
  }
  return(list(loglik, res$probs))
}

# The scalars that describe one ODE segment: k events so far, s open clades, and the bound tri of the stored states;
# recol_dd: whether recolonisation is diversity dependent (see DAISIE_recol_dd()).
DAISIE_IW_direct_pars <- function(pars1, Kprime, k, lxm, lxe, s, tri = Inf) {
  list(lac = pars1[1], mu = pars1[2], Kprime = Kprime, gam = pars1[4], laa = pars1[5], M = pars1[6],
       k = k, lxm = lxm, lxe = lxe, s = s, tri = tri, recol_dd = DAISIE_recol_dd())
}

# Called by DAISIE_loglik_IW() with rhs = "direct" or "compact" (compact = TRUE), after the checks and the settings of the
# dimensions that the versions share.  The steps are those of the loop in DAISIE_loglik_IW().
DAISIE_loglik_IW_direct <- function(pars1,
                                    pars2,
                                    datalist,
                                    brts,
                                    clade,
                                    event,
                                    col,
                                    M,
                                    Kprime,
                                    lxm,
                                    lxe,
                                    cond,
                                    methode,
                                    abstolint,
                                    reltolint,
                                    verbose,
                                    compact = FALSE) {
  nbrts <- length(brts)
  if (compact) {
    # the largest number of stored states over the segments (the event steps give the states of the next segment)
    s_seg <- DAISIE_IW_segments(event, nbrts)
    P_seg <- vapply(seq_along(s_seg) - 1, function(k) DAISIE_IW_layout(lxm, lxe, DAISIE_IW_tri(Kprime, k, TRUE))$P, numeric(1))
    kbig <- which.max(P_seg * 2 ^ s_seg)
    DAISIE_IW_check_dim(lxm, lxe, 2 ^ s_seg[kbig], P = P_seg[kbig])
  } else {
    smax <- DAISIE_IW_max_open_clades(event, nbrts)
    DAISIE_IW_check_dim(lxm, lxe, 2 ^ smax)
  }
  s <- 0                                            # number of open clades; sysdim = 2^s
  P <- DAISIE_IW_layout(lxm, lxe, DAISIE_IW_tri(Kprime, 0, compact))$P      # states per configuration
  probs <- rep(0, P)
  probs[1] <- 1                                     # the state (0, 0) of configuration 0 is the first of the vector
  loglik <- 0
  expandvec <- NULL
  for (k in 0:(nbrts - 2)) {
    if (isTRUE(identical(pars2[4], 3))) {
      message(paste('k = ', k, ', sysdim = ', 2 ^ s, sep = ''))
    }
    tri <- DAISIE_IW_tri(Kprime, k, compact)
    P <- DAISIE_IW_layout(lxm, lxe, tri)$P
    iw_parms <- DAISIE_IW_direct_pars(pars1, Kprime, k, lxm, lxe, s, tri)
    probs <- .Call("daisie_odeint_iw_direct", probs, brts[(k + 1):(k + 2)], iw_parms, methode, abstolint, reltolint)
    cp <- DAISIE_IW_check_probs(probs, loglik, verbose); loglik <- cp[[1]]; probs <- cp[[2]]
    if (k < (nbrts - 2)) {
      if (event[k + 2] == 1) {                      # colonization
        # (the branch of DAISIE_loglik_IW() for a recolonization without speciation can not be reached there:
        # test_for_colonization is TRUE in both cases)
        probs <- .Call("daisie_iw_event", probs, 1L, iw_parms, 0L)
        s <- s + 1
        expandvec <- c(expandvec, clade[k + 2])
      } else if (event[k + 2] == 2) {               # first speciation in clade: its bit is removed
        tocollapse <- which(expandvec == clade[k + 2])
        if (length(tocollapse) != 1) {
          stop("Internal error: the clade of a first speciation is not open.")
        }
        probs <- .Call("daisie_iw_event", probs, 3L, iw_parms, as.integer(tocollapse - 1))
        s <- s - 1
        expandvec <- expandvec[-tocollapse]
      } else {                                      # speciation
        probs <- .Call("daisie_iw_event", probs, 2L, iw_parms, 0L)
      }
      cp <- DAISIE_IW_check_probs(probs, loglik, verbose); loglik <- cp[[1]]; probs <- cp[[2]]
    }
  }
  expandedclades <- which(pracma::histc(clade, 1:length(clade))$cnt == 1)
  lexpandedclades <- length(expandedclades)
  status <- rep(0, lexpandedclades)
  if (lexpandedclades > 0) {
    for (i in lexpandedclades:1) {
      if (datalist[[1 + expandedclades[i]]]$stac == 2) {
        status[i] <- 1
      }
    }
  }
  if (length(status) > 0) {
    decstatus <- bin2dec(rev(status))
  } else {
    decstatus <- 0
  }
  final_index <- 1 + P * decstatus                  # the state (0, 0) of configuration decstatus; P of the last segment
  if (final_index > length(probs)) {
    stop("subscript out of bounds")
  }
  loglik <- loglik + log(probs[final_index])

  if (cond > 0) {
    iw_parms <- DAISIE_IW_direct_pars(pars1, Kprime, 0, lxm, lxe, 0, DAISIE_IW_tri(Kprime, 0, compact))
    probs <- rep(0, DAISIE_IW_layout(lxm, lxe, iw_parms$tri)$P)
    probs[1] <- 1
    probs <- .Call("daisie_odeint_iw_direct", probs, c(min(brts), 0), iw_parms, methode, abstolint, reltolint)
    logcond <- log1p(-probs[1])
    if (logcond == -Inf) {
      message('Parameters lead to probability of extinction of 1. Loglik is set to -Inf')
      loglik <- -Inf
    } else {
      loglik <- loglik - logcond
    }
  }
  print_parameters_and_loglik(pars = pars1,
                              loglik = loglik,
                              verbose = pars2[4],
                              type = 'island_loglik')
  return(as.numeric(loglik))
}
