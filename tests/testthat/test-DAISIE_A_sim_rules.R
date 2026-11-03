# Unit tests of the rules of the DAISIE_A (version 10) simulation, on
# hand-made states, for the special cases of the model description: gap pairs,
# stars and compact clusters, anagenesis, splits of a non-endemic, the copy
# rule, refilled gaps and the protracted option. The mainland is node n + 1.

rule_setup <- function(n, conn, k, gam = 0.5, protracted = FALSE) {
  pars <- DAISIE_A_create_pars(n_islands = n, lac = 0.5, mu = 0.5, gam = gam,
                               alpha = 1, k = k, disp = 1,
                               connectivity = conn)
  cfg <- DAISIE_A_core_cfg(pars, 10, protracted)
  st <- list(lin = list(),
             gen = matrix(NA_real_, 64, 8, dimnames = list(NULL, .DA_GEN_COLS)),
             n_gen = 0L, n_clades = 1L, maxspecID = 1, mainland_n = 1)
  st <- da_add_node(st, parent = 0, time = 1, kind = 1, side = 0, clade = 1,
                    species = 1, time_split = 1)
  list(cfg = cfg, st = st)
}
rule_lineage <- function(cfg, type, islands, pairs = list()) {
  L <- da_new_lineage(id = if (type == 1) 1 else 2, anc = 1, colt = 1,
                      type = type, ana = 0, node = 1, side = 0, clade = 1,
                      n = cfg$n)
  L$occ[islands] <- 1L
  for (p in pairs) {
    L <- da_set_pair(L, p[1], p[2], p[3], if (p[3] > 0) p[4] else -Inf)
  }
  L
}
one_event <- function(type, r, a, b = 0L) {
  list(n = 1L, cs = 1, type = type, r = r, a = a, b = b)
}
star4 <- function() {
  C <- matrix(0, 4, 4)
  C[1, 2:4] <- 1
  C[2:4, 1] <- 1
  C
}

test_that("loss of an intermediate island on a line leaves a gap pair", {
  x <- rule_setup(3, "line", k = 3)
  L <- rule_lineage(x$cfg, 3, 1:3, list(c(1, 2, 1, 2.0), c(2, 3, 0, 0)))
  x$st$lin[[1]] <- L
  x$st <- da_apply(x$st, one_event(2L, 1L, 2L), 1L, 3, x$cfg)
  L <- x$st$lin[[1]]
  expect_length(x$st$lin, 1)
  expect_equal(L$occ, c(1L, 0L, 1L))
  # the gap pair starts at the larger of the two states, with its time
  expect_equal(L$S[1, 3], 1)
  expect_equal(L$Tau[1, 3], 2)
  expect_true(all(L$S[2, ] == -1))
  expect_silent(da_check_invariants(x$st, x$cfg))
})

test_that("losing the hub of a star gives gap pairs between all satellites", {
  x <- rule_setup(4, star4(), k = 2)
  L <- rule_lineage(x$cfg, 3, 1:4,
                    list(c(1, 2, 1, 2.0), c(1, 3, 0, 0), c(1, 4, 1, 3.0)))
  x$st$lin[[1]] <- L
  x$st <- da_apply(x$st, one_event(2L, 1L, 1L), 1L, 4, x$cfg)
  L <- x$st$lin[[1]]
  expect_length(x$st$lin, 1)  # no split: all gap pairs start below k
  expect_equal(c(L$S[2, 3], L$S[2, 4], L$S[3, 4]), c(1, 1, 1))
  # equal states: the later step 0 -> 1 is inherited
  expect_equal(c(L$Tau[2, 3], L$Tau[2, 4], L$Tau[3, 4]), c(2, 3, 3))
  expect_silent(da_check_invariants(x$st, x$cfg))
})

test_that("losing the hub of a compact cluster gives a multifurcation", {
  for (prot in c(FALSE, TRUE)) {
    x <- rule_setup(4, "full", k = 1, protracted = prot)
    L <- rule_lineage(x$cfg, 3, 1:4,
                      list(c(1, 2, 0, 0), c(1, 3, 0, 0), c(1, 4, 0, 0),
                           c(2, 3, 1, 2.0), c(2, 4, 1, 2.5), c(3, 4, 1, 3.0)))
    x$st$lin[[1]] <- L
    set.seed(1)
    x$st <- da_apply(x$st, one_event(2L, 1L, 1L), 1L, 4, x$cfg)
    expect_length(x$st$lin, 3)
    ranges <- sort(vapply(x$st$lin, function(l) which(l$occ == 1L), 0))
    expect_equal(ranges, 2:4)
    expect_true(all(vapply(x$st$lin, function(l) l$type, 0) == 3))
    expect_equal(length(unique(vapply(x$st$lin, function(l) l$id, 0))), 3)
    gen <- x$st$gen[seq_len(x$st$n_gen), , drop = FALSE]
    nodes <- gen[gen[, "kind"] == 2, , drop = FALSE]
    # two binary nodes, the second a child of the first: a zero-length branch
    expect_equal(nrow(nodes), 2)
    expect_equal(unname(nodes[2, "parent"]), unname(nodes[1, "node"]))
    expect_equal(unname(nodes[, "time_split"]), c(4, 4))
    if (prot) {
      # dated at the latest step 0 -> 1 among the cut pairs, never before the
      # parent node
      expect_true(all(nodes[, "time"] <= 4))
      expect_true(nodes[2, "time"] >= nodes[1, "time"])
    } else {
      expect_equal(unname(nodes[, "time"]), c(4, 4))
    }
    expect_silent(da_check_invariants(x$st, x$cfg))
  }
})

test_that("divergence from the mainland gives anagenesis", {
  for (prot in c(FALSE, TRUE)) {
    x <- rule_setup(1, "full", k = 2, protracted = prot)
    x$st$lin[[1]] <- rule_lineage(x$cfg, 1, 1, list(c(1, 2, 1, 1.5)))
    x$st <- da_apply(x$st, one_event(5L, 1L, 1L, 2L), 1L, 2, x$cfg)
    L <- x$st$lin[[1]]
    expect_equal(c(L$type, L$ana, L$id), c(2, 1, 2))
    expect_true(all(L$S == -1))
    # no node in the genealogy; with the protracted option the colonisation is
    # dated at the start of divergence from the mainland
    expect_equal(x$st$n_gen, 1L)
    expect_equal(L$colt, if (prot) 1.5 else 1)
    expect_equal(unname(x$st$gen[1, "time"]), if (prot) 1.5 else 1)
  }
})

test_that("a non-endemic splits: the part with the mainland stays non-endemic", {
  x <- rule_setup(2, "line", k = 1, gam = c(0.5, 0))
  x$st$lin[[1]] <- rule_lineage(x$cfg, 1, 1:2, list(c(1, 3, 0, 0),
                                                    c(1, 2, 0, 0)))
  x$st <- da_apply(x$st, one_event(5L, 1L, 1L, 2L), 1L, 2, x$cfg)
  expect_length(x$st$lin, 2)
  ne <- x$st$lin[[1]]
  en <- x$st$lin[[2]]
  expect_equal(c(ne$type, ne$id), c(1, 1))
  expect_equal(ne$occ, c(1L, 0L))
  expect_equal(ne$S[1, 3], 0)
  expect_equal(c(en$type, en$side), c(3, 2))
  expect_equal(en$occ, c(0L, 1L))
  expect_true(all(en$S == -1))
  expect_equal(unname(x$st$gen[2, c("parent", "kind")]), c(1, 2))
})

test_that("founders copy their source population", {
  x <- rule_setup(3, "line", k = 3)
  x$st$lin[[1]] <- rule_lineage(x$cfg, 1, 1, list(c(1, 4, 1, 1.2)))
  # dispersal 1 -> 2: d12 = 0 and e_2 = e_1
  x$st <- da_apply(x$st, one_event(4L, 1L, 1L, 2L), 1L, 2, x$cfg)
  L <- x$st$lin[[1]]
  expect_equal(c(L$S[1, 2], L$S[2, 4], L$Tau[2, 4]), c(0, 1, 1.2))
  # colonisation of island 3 from the mainland: e_3 = 0 and d23 = e_2
  x$st <- da_apply(x$st, one_event(1L, 1L, 3L), 1L, 3, x$cfg)
  L <- x$st$lin[[1]]
  expect_equal(c(L$S[3, 4], L$S[2, 3], L$Tau[2, 3]), c(0, 1, 1.2))
  expect_equal(L$occ, c(1L, 1L, 1L))
  expect_silent(da_check_invariants(x$st, x$cfg))
  # gene flow resets a pair; recolonisation resets a mainland pair
  x$st <- da_apply(x$st, one_event(4L, 1L, 3L, 2L), 1L, 4, x$cfg)
  expect_equal(x$st$lin[[1]]$S[2, 3], 0)
  x$st <- da_apply(x$st, one_event(1L, 1L, 2L), 1L, 4, x$cfg)
  expect_equal(x$st$lin[[1]]$S[2, 4], 0)
  expect_equal(x$st$lin[[1]]$colt, 1)  # no reset of the colonisation time
})

test_that("a refilled gap is removed when its islands are connected again", {
  x <- rule_setup(3, "line", k = 3)
  x$st$lin[[1]] <- rule_lineage(x$cfg, 3, c(1, 3), list(c(1, 3, 1, 2.0)))
  x$st <- da_apply(x$st, one_event(4L, 1L, 1L, 2L), 1L, 3, x$cfg)
  L <- x$st$lin[[1]]
  expect_equal(c(L$S[1, 2], L$S[2, 3], L$Tau[2, 3]), c(0, 1, 2))
  expect_equal(L$S[1, 3], -1)
  expect_silent(da_check_invariants(x$st, x$cfg))
})

test_that("a gap pair stays while the new route has a pair at k", {
  # non-endemic on islands 3 and 4 of a star with hub 1; the mainland does not
  # reach island 3, which is held in the lineage by its gap pair with island 4
  x <- rule_setup(4, star4(), k = 2, gam = c(0.4, 0.1, 0, 0.1))
  x$st$lin[[1]] <- rule_lineage(x$cfg, 1, 3:4,
                                list(c(3, 4, 0, 0), c(3, 5, 2, 1.1),
                                     c(4, 5, 1, 1.3)))
  # the hub is colonised from the mainland: its pair with island 3 copies
  # (mainland, 3) and is at k, so the gap pairs of island 3 are kept
  x$st <- da_apply(x$st, one_event(1L, 1L, 1L), 1L, 2, x$cfg)
  L <- x$st$lin[[1]]
  expect_length(x$st$lin, 1)
  expect_equal(c(L$S[1, 3], L$S[1, 4], L$S[1, 5]), c(2, 1, 0))
  expect_equal(c(L$S[3, 4], L$S[3, 5]), c(0, 2))
  expect_silent(da_check_invariants(x$st, x$cfg))
})

test_that("protracted within-island speciation of a widespread lineage", {
  x <- rule_setup(2, "full", k = 1, protracted = TRUE)
  x$st$lin[[1]] <- rule_lineage(x$cfg, 1, 1:2, list(c(1, 2, 0, 0),
                                                    c(1, 3, 0, 0),
                                                    c(2, 3, 0, 0)))
  x$st <- da_apply(x$st, one_event(3L, 1L, 1L), 1L, 2, x$cfg)
  expect_length(x$st$lin, 3)
  expect_equal(x$st$lin[[1]]$occ, c(0L, 1L))
  expect_equal(x$st$lin[[1]]$type, 1)
  expect_true(all(vapply(x$st$lin[2:3], function(l) l$occ[1] == 1L &&
                           l$type == 3, TRUE)))
  gen <- x$st$gen[seq_len(x$st$n_gen), , drop = FALSE]
  # node 1 (island-1 population vs the rest) and node 2 (the two new species)
  # coincide for k = 1: a zero-length branch
  expect_equal(nrow(gen), 3)
  expect_equal(unname(gen[3, "parent"]), 2)
  expect_equal(unname(gen[2:3, "time"]), c(2, 2))
  # for k = 2 the within-island state must step twice; gene flow resets it
  x <- rule_setup(2, "full", k = 2, protracted = TRUE)
  x$st$lin[[1]] <- rule_lineage(x$cfg, 3, 1:2, list(c(1, 2, 0, 0)))
  x$st <- da_apply(x$st, one_event(3L, 1L, 1L), 1L, 2, x$cfg)
  expect_equal(c(x$st$lin[[1]]$s[1], x$st$lin[[1]]$stau[1]), c(1, 2))
  ev <- da_events(x$st, x$cfg, 2)
  expect_true(any(ev$type == 4 & ev$a == 2 & ev$b == 1))  # gene flow into 1
  x$st <- da_apply(x$st, one_event(4L, 1L, 2L, 1L), 1L, 2.5, x$cfg)
  expect_equal(x$st$lin[[1]]$s[1], 0)
  x$st <- da_apply(x$st, one_event(3L, 1L, 1L), 1L, 3, x$cfg)
  x$st <- da_apply(x$st, one_event(3L, 1L, 1L), 1L, 4, x$cfg)
  gen <- x$st$gen[seq_len(x$st$n_gen), , drop = FALSE]
  # node 1 dated at the step 0 -> 1 (time 3), node 2 at the completion (4)
  expect_equal(unname(gen[2:3, "time"]), c(3, 4))
  expect_equal(unname(gen[2:3, "time_split"]), c(4, 4))
})

test_that("only events that add a species are diversity dependent", {
  pars <- DAISIE_A_create_pars(n_islands = 2, lac = 0.5, mu = 0, K = 2,
                               gam = 0.3, alpha = 1, k = 2, disp = 0.7)
  cfg <- DAISIE_A_core_cfg(pars, 10, FALSE)
  st <- list(lin = list(), gen = matrix(NA_real_, 8, 8), n_gen = 0L,
             n_clades = 1L, maxspecID = 2, mainland_n = 1)
  # the non-endemic on both islands (pairs at 1) and an endemic on island 1
  L <- rule_lineage(cfg, 1, 1:2, list(c(1, 2, 1, 1), c(1, 3, 1, 1),
                                      c(2, 3, 1, 1)))
  E <- rule_lineage(cfg, 3, 1)
  st$lin <- list(L, E)
  ev <- da_events(st, cfg, 1)
  rate <- function(type, r, a, b = 0L) {
    i <- which(ev$type == type & ev$r == r & ev$a == a & ev$b == b)
    if (length(i) == 0) return(0)
    ev$cs[i] - c(0, ev$cs)[i]
  }
  f <- 1 - c(2, 1) / 2
  expect_equal(rate(1L, 1L, 1L), 0.3)              # recolonisation: not DD
  expect_equal(rate(3L, 1L, 2L), 0.5 * f[2])       # cladogenesis on 2
  expect_equal(rate(3L, 1L, 1L), 0)                # island 1 is full
  expect_equal(rate(4L, 2L, 1L, 2L), 0.7 * f[2])   # range expansion
  expect_equal(rate(4L, 1L, 1L, 2L), 0.7)          # gene flow: not DD
  expect_equal(sum(ev$type == 5), 3)               # three pairs below k
})
