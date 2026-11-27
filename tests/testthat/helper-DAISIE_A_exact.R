# Exact distribution of the state of the DAISIE_A model (version 10,
# divergence states, protracted = FALSE) for one mainland species on a small
# archipelago, from the forward (master) equation.
#
# This is an independent implementation of the model rules (written from the
# model description, not from the simulation code): the states reachable from
# the empty archipelago are enumerated, the generator matrix of the
# continuous-time Markov chain is built, and the distribution at time T is
# computed by uniformization. The diversity-dependence factors bound the
# numbers of species per island, so the state space is finite.
#
# A state is a list of lineages. A lineage is a list with
#   type   1 non-endemic (I), 2 anagenetic endemic (A), 3 cladogenetic (C);
#   occ    0/1 vector of length n (its range);
#   S      (n + 1) x (n + 1) matrix of divergence states of its tracked pairs
#          (node n + 1 is the mainland), -1 = not tracked;
#   clade  colonization event (only the grouping matters).

exact_model <- function(pars) {
  n <- pars$n_islands
  n1 <- n + 1L
  adj <- matrix(FALSE, n1, n1)
  adj[1:n, 1:n] <- pars$adjacency
  adj[n1, 1:n] <- pars$gam > 0
  adj[1:n, n1] <- pars$gam > 0
  diag(adj) <- FALSE
  list(n = n, n1 = n1, adj = adj, lac = pars$lac, mu = pars$mu, K = pars$K,
       gam = pars$gam, alpha = pars$alpha, k = pars$k, D = pars$D,
       split_single = pars$clado_single == "split")
}

exact_nodes <- function(L, m) {
  c(which(L$occ == 1), if (L$type == 1) m$n1)
}

# connected components of the nodes nd; edge(u, w) says whether u and w are
# linked; returns a vector over all nodes (0 for nodes not in nd)
exact_components <- function(nd, edge, m) {
  comp <- integer(m$n1)
  nc <- 0L
  for (v in nd) {
    if (comp[v] != 0L) next
    nc <- nc + 1L
    comp[v] <- nc
    queue <- v
    while (length(queue) > 0) {
      u <- queue[1]
      queue <- queue[-1]
      for (w in nd) {
        if (comp[w] == 0L && edge(u, w)) {
          comp[w] <- nc
          queue <- c(queue, w)
        }
      }
    }
  }
  comp
}

exact_set_pair <- function(L, a, b, s) {
  L$S[a, b] <- s
  L$S[b, a] <- s
  L
}

# drop gap pairs (tracked pairs of non-adjacent nodes) whose nodes are
# connected through adjacent pairs below k
exact_drop_gaps <- function(L, m) {
  nd <- exact_nodes(L, m)
  comp <- exact_components(nd, function(u, w) {
    m$adj[u, w] && L$S[u, w] >= 0 && L$S[u, w] < m$k
  }, m)
  for (a in nd) for (b in nd) {
    if (a < b && L$S[a, b] >= 0 && !m$adj[a, b] && comp[a] == comp[b]) {
      L <- exact_set_pair(L, a, b, -1L)
    }
  }
  L
}

# a new population on j founded from node src: the copy rule (pairs of j with
# the adjacent nodes l of the lineage start at the state of src with l, or at
# the minimax path value if src and l have no tracked pair)
exact_add_island <- function(L, j, src, m) {
  nd <- exact_nodes(L, m)
  d <- matrix(Inf, m$n1, m$n1)
  for (a in nd) for (b in nd) {
    if (a == b) d[a, b] <- 0 else if (L$S[a, b] >= 0) d[a, b] <- L$S[a, b]
  }
  for (v in nd) for (a in nd) for (b in nd) {
    d[a, b] <- min(d[a, b], max(d[a, v], d[v, b]))
  }
  for (l in nd) {
    if (!m$adj[j, l]) next
    s <- if (l == src) 0L else if (L$S[src, l] >= 0) L$S[src, l] else d[src, l]
    L <- exact_set_pair(L, j, l, as.integer(s))
  }
  L$occ[j] <- 1L
  exact_drop_gaps(L, m)
}

# loss of the population on l: its pairs disappear; two neighbors of l that
# are no longer connected by the remaining pairs get a gap pair at the larger
# of their states with l
exact_remove_island <- function(L, l, m) {
  nd <- exact_nodes(L, m)
  nb <- nd[nd != l & L$S[l, nd] >= 0]
  ls <- L$S[l, ]
  L$S[l, ] <- -1L
  L$S[, l] <- -1L
  L$occ[l] <- 0L
  nd2 <- setdiff(nd, l)
  comp <- exact_components(nd2, function(u, w) L$S[u, w] >= 0, m)
  if (length(nb) >= 2) {
    for (x in 1:(length(nb) - 1)) for (y in (x + 1):length(nb)) {
      a <- nb[x]
      b <- nb[y]
      if (comp[a] != comp[b]) L <- exact_set_pair(L, a, b, max(ls[a], ls[b]))
    }
  }
  L
}

exact_restrict <- function(L, P, m) {
  out <- setdiff(seq_len(m$n1), P)
  L$S[out, ] <- -1L
  L$S[, out] <- -1L
  L$occ[setdiff(seq_len(m$n), P)] <- 0L
  L
}

# the split rule: the lineages that replace L (a list with L itself if the
# graph of its pairs below k is connected)
exact_split <- function(L, m) {
  nd <- exact_nodes(L, m)
  comp <- exact_components(nd, function(u, w) {
    L$S[u, w] >= 0 && L$S[u, w] < m$k
  }, m)
  pieces <- split(nd, comp[nd])
  if (length(pieces) == 1) return(list(L))
  new_endemic <- function(P) {
    B <- exact_restrict(L, P, m)
    B$type <- 3L
    B
  }
  if (L$type == 1) {
    has_m <- vapply(pieces, function(P) m$n1 %in% P, logical(1))
    main <- pieces[[which(has_m)]]
    others <- pieces[!has_m]
    if (length(main) > 1) {
      return(c(list(exact_restrict(L, main, m)), lapply(others, new_endemic)))
    }
    # the mainland is alone: anagenesis with one other piece, otherwise one
    # new endemic per piece
    if (length(others) == 1) {
      A <- exact_restrict(L, others[[1]], m)
      A$type <- 2L
      return(list(A))
    }
    return(lapply(others, new_endemic))
  }
  lapply(pieces, new_endemic)
}

exact_new_lineage <- function(type, i, clade, m) {
  L <- list(type = as.integer(type), occ = integer(m$n),
            S = matrix(-1L, m$n1, m$n1), clade = clade)
  L$occ[i] <- 1L
  L
}

# a cladogenetic species that is the only species left of its colonization
# event becomes anagenetic ("Clado_extinct")
exact_clade_singleton <- function(lins, clade) {
  in_clade <- which(vapply(lins, function(L) L$clade == clade, logical(1)))
  if (length(in_clade) == 1 && lins[[in_clade]]$type == 3) {
    lins[[in_clade]]$type <- 2L
  }
  lins
}

exact_lineage_key <- function(L) {
  paste0(L$type, "/", paste(L$occ, collapse = ""), "/",
         paste(L$S[upper.tri(L$S)], collapse = ","))
}

# canonical form: clades relabeled and lineages sorted
exact_canonical <- function(lins) {
  if (length(lins) == 0) return(list(key = "<empty>", lins = list()))
  keys <- vapply(lins, exact_lineage_key, "")
  clades <- vapply(lins, function(L) L$clade, numeric(1))
  ck <- tapply(keys, clades, function(x) paste(sort(x), collapse = "+"))
  ord_clades <- names(ck)[order(ck)]
  out <- list()
  for (q in seq_along(ord_clades)) {
    idx <- which(clades == as.numeric(ord_clades[q]))
    idx <- idx[order(keys[idx])]
    for (r in idx) {
      L <- lins[[r]]
      L$clade <- q
      out[[length(out) + 1]] <- L
    }
  }
  list(key = paste(sort(ck), collapse = "|"), lins = out)
}

# all transitions out of a state: a list of (rate, new list of lineages)
exact_transitions <- function(lins, m) {
  n <- m$n
  out <- list()
  add <- function(rate, new_lins) {
    if (rate > 0) out[[length(out) + 1]] <<- list(rate = rate, lins = new_lins)
  }
  occ_count <- numeric(n)
  for (L in lins) occ_count <- occ_count + L$occ
  f <- pmax(0, 1 - occ_count / m$K)
  new_clade <- if (length(lins) == 0) 1 else
    max(vapply(lins, function(L) L$clade, numeric(1))) + 1
  ne <- which(vapply(lins, function(L) L$type == 1, logical(1)))
  # colonization from the mainland
  for (i in seq_len(n)) {
    if (m$gam[i] == 0) next
    if (length(ne) == 0) {
      L <- exact_new_lineage(1, i, new_clade, m)
      L <- exact_set_pair(L, i, m$n1, 0L)
      add(m$gam[i] * f[i], c(lins, list(L)))
      next
    }
    L <- lins[[ne]]
    if (L$occ[i] == 0) {
      new_lins <- lins
      new_lins[[ne]] <- exact_add_island(L, i, m$n1, m)
      add(m$gam[i] * f[i], new_lins)
    } else if (sum(L$occ) == 1) {
      # recolonization of the only island: a new colonization event
      N <- exact_new_lineage(1, i, new_clade, m)
      N <- exact_set_pair(N, i, m$n1, 0L)
      new_lins <- c(lins[-ne], list(N))
      new_lins <- exact_clade_singleton(new_lins, L$clade)
      add(m$gam[i], new_lins)
    } else if (L$S[i, m$n1] > 0) {
      new_lins <- lins
      new_lins[[ne]] <- exact_drop_gaps(exact_set_pair(L, i, m$n1, 0L), m)
      add(m$gam[i], new_lins)
    }
  }
  for (r in seq_along(lins)) {
    L <- lins[[r]]
    R <- which(L$occ == 1)
    # local extinction
    for (i in R) {
      if (m$mu[i] == 0) next
      if (length(R) == 1) {
        new_lins <- exact_clade_singleton(lins[-r], L$clade)
      } else {
        new_lins <- c(lins[-r], exact_split(exact_remove_island(L, i, m), m))
      }
      add(m$mu[i], new_lins)
    }
    # within-island speciation
    for (i in R) {
      rate <- m$lac[i] * f[i]
      if (rate == 0) next
      if (length(R) == 1 && m$split_single) {
        new_lins <- c(lins[-r], list(exact_new_lineage(3, i, L$clade, m),
                                     exact_new_lineage(3, i, L$clade, m)))
      } else {
        new_lins <- c(lins, list(exact_new_lineage(3, i, L$clade, m)))
      }
      add(rate, new_lins)
    }
    # dispersal
    for (i in R) {
      for (j in seq_len(n)) {
        if (j == i || m$D[i, j] == 0) next
        if (L$occ[j] == 0) {
          new_lins <- lins
          new_lins[[r]] <- exact_add_island(L, j, i, m)
          add(m$D[i, j] * f[j], new_lins)
        } else if (L$S[i, j] > 0) {
          new_lins <- lins
          new_lins[[r]] <- exact_drop_gaps(exact_set_pair(L, i, j, 0L), m)
          add(m$D[i, j], new_lins)
        }
      }
    }
    # divergence steps
    nd <- exact_nodes(L, m)
    for (a in nd) for (b in nd) {
      if (a < b && L$S[a, b] >= 0 && L$S[a, b] < m$k) {
        L2 <- exact_set_pair(L, a, b, L$S[a, b] + 1L)
        repl <- if (L2$S[a, b] == m$k) exact_split(L2, m) else list(L2)
        add(m$alpha, c(lins[-r], repl))
      }
    }
  }
  out
}

# enumerate the states reachable from the empty archipelago and build the
# generator matrix
exact_chain <- function(pars, max_states = 50000) {
  m <- exact_model(pars)
  keys <- "<empty>"
  states <- list(list())
  index <- new.env(hash = TRUE)
  index[["<empty>"]] <- 1L
  from <- integer(0)
  to <- integer(0)
  rate <- numeric(0)
  q <- 1L
  while (q <= length(states)) {
    for (tr in exact_transitions(states[[q]], m)) {
      cn <- exact_canonical(tr$lins)
      j <- index[[cn$key]]
      if (is.null(j)) {
        j <- length(states) + 1L
        if (j > max_states) stop("too many states")
        states[[j]] <- cn$lins
        keys[j] <- cn$key
        index[[cn$key]] <- j
      }
      if (j != q) {
        from <- c(from, q)
        to <- c(to, j)
        rate <- c(rate, tr$rate)
      }
    }
    q <- q + 1L
  }
  ns <- length(states)
  Q <- Matrix::sparseMatrix(i = from, j = to, x = rate, dims = c(ns, ns))
  Matrix::diag(Q) <- -Matrix::rowSums(Q)
  list(states = states, keys = keys, Q = Q, model = m)
}

# distribution at time t, starting from the empty archipelago
# (uniformization)
exact_distribution <- function(chain, t, tol = 1e-13) {
  Q <- chain$Q
  ns <- nrow(Q)
  lambda <- max(-Matrix::diag(Q))
  P <- Matrix::Diagonal(ns) + Q / lambda
  v <- c(1, rep(0, ns - 1))
  out <- numeric(ns)
  w <- stats::dpois(0, lambda * t)
  cum <- w
  out <- out + w * v
  j <- 0
  while (cum < 1 - tol) {
    j <- j + 1
    v <- as.numeric(v %*% P)
    w <- stats::dpois(j, lambda * t)
    cum <- cum + w
    out <- out + w * v
  }
  out
}

# observable summary of a state: the species with their type and range,
# grouped by colonization event
exact_summary <- function(lins) {
  if (length(lins) == 0) return("")
  sp <- vapply(lins, function(L) {
    paste0(c("I", "A", "C")[L$type], paste(L$occ, collapse = ""))
  }, "")
  clade <- vapply(lins, function(L) L$clade, numeric(1))
  paste(sort(tapply(sp, clade, function(x) paste(sort(x), collapse = "+"))),
        collapse = ",")
}

# the same summary of a simulated mainland species (output of
# DAISIE_A_call_cpp(): columns 4 type, 8 clade, then the range)
exact_sim_summary <- function(core, n) {
  spec <- core$spec
  if (nrow(spec) == 0) return("")
  occ <- spec[, 8 + seq_len(n), drop = FALSE]
  sp <- paste0(c("I", "A", "C")[spec[, 4]], apply(occ, 1, paste, collapse = ""))
  paste(sort(tapply(sp, spec[, 8], function(x) paste(sort(x), collapse = "+"))),
        collapse = ",")
}
