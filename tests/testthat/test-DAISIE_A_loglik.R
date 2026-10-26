tip_clade <- function(col_time, range, type = "E") {
  list(colonisation_time = col_time,
       nodes = data.frame(node = 1, parent = 0, age = 0, range = range,
                          type = type))
}

random_clade <- function(col_time, node_ages, range = "1") {
  node_ages <- sort(node_ages, decreasing = TRUE)
  nodes <- data.frame(node = 1, parent = 0, age = node_ages[1], range = NA,
                      type = NA)
  open <- c(1, 1)
  id <- 1
  for (a in node_ages[-1]) {
    k <- sample.int(length(open), 1)
    id <- id + 1
    nodes <- rbind(nodes, data.frame(node = id, parent = open[k], age = a,
                                     range = NA, type = NA))
    open <- c(open[-k], id, id)
  }
  for (p in open) {
    id <- id + 1
    nodes <- rbind(nodes, data.frame(node = id, parent = p, age = 0,
                                     range = range, type = "E"))
  }
  list(colonisation_time = col_time, nodes = nodes)
}

test_that("with one island the likelihood equals DAISIE_DE", {
  p <- c(0.8, 0.4, 0.4, 0.06, 0.5)
  pars <- DAISIE_A_create_pars(1, lac = p[1], mu = p[2], gam = p[4],
                               laa = p[5])
  t0 <- 5
  # island-level anagenesis does not apply to a single island
  pars_ia <- DAISIE_A_create_pars(1, lac = p[1], mu = p[2], gam = p[4],
                                  laa = p[5], laa_isl = 0.3)
  el <- list(clades = list(random_clade(3.1, c(1.9, 1.2, 0.7))),
             nonendemic_range = "1")
  expect_equal(DAISIE_A_loglik_element(pars_ia, el, t0),
               DAISIE_A_loglik_element(pars, el, t0))
  de <- function(brts, stac, miss, sampling) {
    DAISIE_DE_loglik(pars1 = p, brts = brts, missnumspec = miss, stac = stac,
                     methode = "odeint::runge_kutta_cash_karp54",
                     reltolint = 1e-14, abstolint = 1e-14,
                     sampling = if (sampling == "exact") "n" else "rho")
  }
  expect_equal(
    DAISIE_A_loglik_element(pars, NULL, t0),
    DAISIE_DE_logp0(t0, p, methode = "odeint::runge_kutta_cash_karp54",
                    reltolint = 1e-14, abstolint = 1e-14),
    tolerance = 1e-8)
  set.seed(1)
  for (sampling in c("exact", "rho")) {
    expect_equal(
      DAISIE_A_loglik_element(
        pars, list(clades = list(tip_clade(1.7, "1", "I"))), t0,
        sampling = sampling),
      de(c(t0, 1.7), 4, 0, sampling), tolerance = 1e-8)
    for (m in c(0, 2)) {
      el <- list(clades = list(tip_clade(2.3, "1")), missing = rep("1", m))
      expect_equal(DAISIE_A_loglik_element(pars, el, t0, sampling = sampling),
                   de(c(t0, 2.3), 2, m, sampling), tolerance = 1e-8)
      ages <- c(1.9, 1.2, 0.7)
      el <- list(clades = list(random_clade(3.1, ages)), missing = rep("1", m),
                 nonendemic_range = "1")
      expect_equal(DAISIE_A_loglik_element(pars, el, t0, sampling = sampling),
                   de(c(t0, 3.1, ages), 3, m, sampling), tolerance = 1e-8)
    }
  }
})

test_that("internal consistency of the likelihood", {
  pars <- DAISIE_A_create_pars(2, lac = c(0.4, 0.6), mu = c(0.3, 0.2),
                               gam = c(0.06, 0.04), laa = 0.5, lav = 0.4,
                               disp = 0.3, laa_isl = c(0.2, 0.1))
  t0 <- 1.5
  # a missing species of unknown range = sum over its possible ranges
  el <- list(clades = list(tip_clade(1.2, "10")), missing = NA)
  known <- vapply(c("10", "01", "11"), function(r) {
    el$missing <- r
    exp(DAISIE_A_loglik_element(pars, el, t0))
  }, numeric(1))
  el$missing <- NA
  expect_equal(exp(DAISIE_A_loglik_element(pars, el, t0)), sum(known),
               tolerance = 1e-7)
  # rho sampling with a complete phylogeny equals exact sampling
  el2 <- list(clades = list(tip_clade(1.2, "11", "I")))
  expect_equal(DAISIE_A_loglik_element(pars, el2, t0, sampling = "rho"),
               DAISIE_A_loglik_element(pars, el2, t0, sampling = "exact"))
  # a non-endemic with unknown colonisation time integrates over the
  # colonisation time
  f <- function(tc) {
    vapply(tc, function(t) {
      exp(DAISIE_A_loglik_element(
        pars, list(clades = list(tip_clade(t, "11", "I"))), t0))
    }, numeric(1))
  }
  expect_equal(
    exp(DAISIE_A_loglik_element(pars, list(nonendemic_range = "11"), t0)),
    stats::integrate(f, 0, t0, rel.tol = 1e-7)$value, tolerance = 1e-6)
  # the likelihood does not depend on the order of the children
  nodes <- data.frame(node = 1:3, parent = c(0, 1, 1), age = c(0.6, 0, 0),
                      range = c(NA, "10", "11"), type = c(NA, "E", "I"))
  nodes2 <- nodes
  nodes2[2:3, c("range", "type")] <- nodes[3:2, c("range", "type")]
  expect_equal(
    DAISIE_A_loglik_element(pars, list(clades = list(
      list(colonisation_time = 1.2, nodes = nodes))), t0),
    DAISIE_A_loglik_element(pars, list(clades = list(
      list(colonisation_time = 1.2, nodes = nodes2))), t0))
})

test_that("simulated data can be converted and evaluated", {
  pars <- DAISIE_A_create_pars(3, lac = 0.6, mu = 0.3, gam = 0.05, laa = 0.4,
                               lav = 0.3, disp = 0.4, connectivity = "line")
  set.seed(4)
  sim <- DAISIE_A_sim(time = 2, M = 15, pars = pars, replicates = 1,
                      verbose = FALSE)
  dl <- DAISIE_A_sim_to_datalist(sim[[1]])
  expect_equal(length(dl), length(sim[[1]]))
  ll <- DAISIE_A_loglik(pars, dl)
  expect_true(is.finite(ll))
  expect_true(DAISIE_A_loglik(pars, dl, cond = 1) > ll)
  set.seed(5)
  dl_miss <- DAISIE_A_sim_to_datalist(sim[[1]], prop_missing = 0.5)
  expect_true(is.finite(DAISIE_A_loglik(pars, dl_miss)))
  expect_true(is.finite(DAISIE_A_loglik(pars, dl_miss, sampling = "rho")))
})

test_that("pruning and conversion of phylogenies", {
  nodes <- data.frame(node = 1:5, parent = c(0, 1, 1, 2, 2),
                      age = c(2, 1, 0, 0, 0), range = c(NA, NA, "1", "1", "1"),
                      type = c(NA, NA, "E", "E", "E"))
  pr <- DAISIE_A_prune_nodes(nodes, keep = c(3, 5))
  expect_equal(nrow(pr), 3)
  expect_equal(pr$age, c(2, 0, 0))
  expect_equal(pr$parent, c(0, 1, 1))
  pr2 <- DAISIE_A_prune_nodes(nodes, keep = c(4, 5))
  expect_equal(pr2$age, c(1, 0, 0))
  skip_if_not_installed("ape")
  phy <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
  cl <- DAISIE_A_phylo_to_clade(phy, c(a = "10", b = "11", c = "01"),
                                c(a = "E", b = "E", c = "I"),
                                colonisation_time = 3)
  expect_equal(cl$nodes$age[cl$nodes$parent == 0], 2)
  expect_equal(sort(cl$nodes$range[!is.na(cl$nodes$range)]),
               c("01", "10", "11"))
  pars <- DAISIE_A_create_pars(2, lac = 0.5, mu = 0.2, gam = 0.1, laa = 0.3,
                               lav = 0.2, disp = 0.3)
  expect_true(is.finite(DAISIE_A_loglik_element(
    pars, list(clades = list(cl)), island_age = 4)))
})

test_that("invalid data are rejected", {
  pars <- DAISIE_A_create_pars(2, lac = 0.5, mu = 0.2, gam = 0.1, laa = 0.3)
  bad <- list(clades = list(tip_clade(1, "1")))
  expect_error(DAISIE_A_loglik_element(pars, bad, 2), "range")
  two_I <- list(clades = list(tip_clade(1.5, "10", "I"),
                              tip_clade(1, "01", "I")))
  expect_error(DAISIE_A_loglik_element(pars, two_I, 2), "youngest")
  both <- list(clades = list(tip_clade(1, "10", "I")),
               nonendemic_range = "01")
  expect_error(DAISIE_A_loglik_element(pars, both, 2), "not both")
})

test_that("islands that have not emerged cannot be colonised", {
  pars <- DAISIE_A_create_pars(2, lac = 0.4, mu = 0.3, gam = 0.1, laa = 0.3,
                               lav = 0.2, disp = 0.3, island_ages = c(3, 1))
  t0 <- 3
  # colonisation of island 2 before it emerged is impossible (without
  # dispersal and local extinction, a non-endemic on island 2 alone must have
  # colonised island 2 directly)
  pars0 <- DAISIE_A_create_pars(2, lac = 0.4, mu = 0, gam = 0.1, laa = 0.3,
                                lav = 0.2, disp = 0, island_ages = c(3, 1))
  expect_equal(DAISIE_A_loglik_element(
    pars0, list(clades = list(tip_clade(2, "01", "I"))), t0), -Inf)
  # with dispersal it is possible: colonise island 1, expand, lose island 1
  expect_true(is.finite(DAISIE_A_loglik_element(
    pars, list(clades = list(tip_clade(2, "01", "I"))), t0)))
  expect_true(is.finite(DAISIE_A_loglik_element(
    pars, list(clades = list(tip_clade(0.5, "01", "I"))), t0)))
  # with equal island ages the result is the same as without island ages
  pars_eq <- DAISIE_A_create_pars(2, lac = 0.4, mu = 0.3, gam = 0.1,
                                  laa = 0.3, lav = 0.2, disp = 0.3,
                                  island_ages = c(3, 3))
  pars_no <- DAISIE_A_create_pars(2, lac = 0.4, mu = 0.3, gam = 0.1,
                                  laa = 0.3, lav = 0.2, disp = 0.3)
  el <- list(clades = list(tip_clade(2, "11", "I")))
  expect_equal(DAISIE_A_loglik_element(pars_eq, el, t0),
               DAISIE_A_loglik_element(pars_no, el, t0))
  # the non-endemic with unknown colonisation time integrates over time,
  # also with staggered island ages
  f <- function(tc) {
    vapply(tc, function(t) {
      exp(DAISIE_A_loglik_element(
        pars, list(clades = list(tip_clade(t, "11", "I"))), t0))
    }, numeric(1))
  }
  expect_equal(
    exp(DAISIE_A_loglik_element(pars, list(nonendemic_range = "11"), t0)),
    stats::integrate(f, 0, 1, rel.tol = 1e-8)$value +
      stats::integrate(f, 1, t0, rel.tol = 1e-8)$value, tolerance = 1e-6)
})
