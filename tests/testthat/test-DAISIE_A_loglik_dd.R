# Approximate (mean-field) likelihood for finite carrying capacities

ladder_clade <- function(col_time, node_ages, range = "1") {
  node_ages <- sort(node_ages, decreasing = TRUE)
  nodes <- data.frame(node = 1, parent = 0, age = node_ages[1], range = NA,
                      type = NA)
  open <- c(1, 1)
  id <- 1
  for (a in node_ages[-1]) {
    id <- id + 1
    nodes <- rbind(nodes, data.frame(node = id, parent = open[1], age = a,
                                     range = NA, type = NA))
    open <- c(open[-1], id, id)
  }
  for (p in open) {
    id <- id + 1
    nodes <- rbind(nodes, data.frame(node = id, parent = p, age = 0,
                                     range = range, type = "E"))
  }
  list(colonisation_time = col_time, nodes = nodes)
}

test_that("with infinite or very large K the likelihood is the K = Inf one", {
  pars <- DAISIE_A_create_pars(3, lac = c(0.5, 0.4, 0.6),
                               mu = c(0.2, 0.3, 0.25),
                               gam = c(0.03, 0.02, 0.025), laa = 0.4,
                               lav = 0.2, disp = 0.3, connectivity = "line",
                               laa_isl = 0.1)
  set.seed(3)
  sim <- DAISIE_A_sim(time = 3, M = 60, pars = pars, replicates = 1,
                      verbose = FALSE)
  dl <- DAISIE_A_sim_to_datalist(sim[[1]])
  l_inf <- DAISIE_A_loglik(pars, dl)
  big <- pars
  big$K <- rep(1e9, 3)
  # the rate components add up to the full rates
  expect_equal(DAISIE_A_loglik(big, dl), l_inf, tolerance = 1e-8)
  # finite K can be ignored
  small <- pars
  small$K <- c(8, 6, 10)
  expect_equal(DAISIE_A_loglik(small, dl, dd_method = "ignore"), l_inf)
  l_mf <- DAISIE_A_loglik(small, dl)
  expect_true(is.finite(l_mf))
  expect_false(isTRUE(all.equal(l_mf, l_inf)))
  # C++ and deSolve integration agree
  skip_on_cran()
  expect_equal(DAISIE_A_loglik(small, dl, methode = "lsoda"), l_mf,
               tolerance = 1e-6)
})

test_that("for one island the mean-field likelihood is close to DAISIE", {
  skip_on_cran()
  t0 <- 5
  ages <- c(1.9, 1.2, 0.7)
  cases <- list(
    list(el = list(clades = list(ladder_clade(3.1, ages))),
         brts = c(t0, 3.1, ages), stac = 2),
    list(el = list(clades = list(list(
      colonisation_time = 2.3,
      nodes = data.frame(node = 1, parent = 0, age = 0, range = "1",
                         type = "E")))),
      brts = c(t0, 2.3), stac = 2),
    list(el = NULL, brts = c(t0, 0), stac = 0))
  for (K in c(20, 10)) {
    pars <- DAISIE_A_create_pars(1, lac = 0.8, mu = 0.4, gam = 0.06,
                                 laa = 0.5, K = K)
    for (cs in cases) {
      exact <- DAISIE:::DAISIE_loglik(
        pars1 = c(0.8, 0.4, K, 0.06, 0.5), pars2 = c(100, 11, 0, 0),
        brts = cs$brts, stac = cs$stac, missnumspec = 0,
        methode = "odeint::runge_kutta_cash_karp54",
        CS_version = list(model = 1, function_to_optimize = "DAISIE"),
        abstolint = 1e-16, reltolint = 1e-12)
      mf <- DAISIE_A_loglik_element(pars, cs$el, t0)
      ign <- DAISIE_A_loglik_element(pars, cs$el, t0, dd_method = "ignore")
      expect_lt(abs(mf - exact), 0.025)
      # and much closer than ignoring K
      expect_lt(abs(mf - exact), abs(ign - exact) / 3)
    }
  }
})

test_that("more species than K on an island have probability zero", {
  pars <- DAISIE_A_create_pars(2, lac = 0.8, mu = 0.4, gam = 0.06,
                               laa = 0.5, disp = 0.2, K = c(3, 10))
  el <- list(clades = list(ladder_clade(3.1, c(1.9, 1.2, 0.7), "10")))
  expect_equal(DAISIE_A_loglik_element(pars, el, 5), -Inf)
  el$clades[[1]]$nodes$range[!is.na(el$clades[[1]]$nodes$range)] <- "01"
  expect_true(is.finite(DAISIE_A_loglik_element(pars, el, 5)))
  out <- DAISIE_A_loglik_element(pars, el, 5, return_dd = TRUE)
  expect_true(out$dd$converged)
  # reconstructed lineages are counted on the islands of their descendants
  expect_equal(DAISIE:::DAISIE_A_dd_l(out$dd, 0.1), c(0, 4))
  expect_equal(DAISIE:::DAISIE_A_dd_l(out$dd, 3), c(0, 1))
})

test_that("diversity dependence reduces the rates as in the simulation", {
  pars <- DAISIE_A_create_pars(2, lac = c(0.5, 0.7), mu = 0.3,
                               gam = c(0.05, 0.02), laa = 0.4, disp = 0.2,
                               lav = 0.3, K = c(5, 7))
  comps <- DAISIE:::DAISIE_A_dd_components(pars)
  expect_length(comps, 5)
  # the components add up to the rates
  expect_equal(Reduce(`+`, lapply(comps, function(p) p$lac)), pars$lac)
  expect_equal(Reduce(`+`, lapply(comps, function(p) p$gam)), pars$gam)
  expect_equal(Reduce(`+`, lapply(comps, function(p) p$D)), pars$D)
  expect_equal(Reduce(`+`, lapply(comps, function(p) p$mu)), pars$mu)
  expect_equal(Reduce(`+`, lapply(comps, function(p) p$lav)), pars$lav)
  dd <- list(K = pars$K, gt = c(0, 1), gh = matrix(c(1, 2, 1, 2), 2),
             gp = matrix(0, 2, 2))
  W <- DAISIE:::DAISIE_A_dd_wmat(dd, 0.5, l = c(2, 3), own = c(1, 0))
  # nu = (3, 5): a lineage counts itself when it speciates on its own island
  # and not when it disperses to another island; a reconstructed lineage is
  # counted without its own contribution (own) to l
  expect_equal(unname(W[2, "E"]), 1 - 4 / 5)
  expect_equal(unname(W[2, "D"]), 1 - 3 / 5)
  expect_equal(unname(W[4, "E"]), 1 - 3 / 5)
  expect_equal(unname(W[4, "D"]), 1 - 2 / 5)
  expect_equal(unname(W[5, "A"]), 1 - 5 / 7)
  expect_equal(attr(W, "reset"), c(1 - 3 / 5, 1 - 6 / 7))
  expect_equal(unname(W[1, ]), rep(1, 4))
})

test_that("a lineage does not block its own dispersal", {
  # a widespread non-endemic species with K = 1 on both islands: it must
  # have dispersed (or colonised) when it was alone on the archipelago
  pars <- DAISIE_A_create_pars(2, lac = 0.5, mu = 0.2, gam = 0.05,
                               laa = 0.3, disp = 0.4, K = c(1, 1))
  el <- list(clades = list(list(
    colonisation_time = 1.5,
    nodes = data.frame(node = 1, parent = 0, age = 0, range = "11",
                       type = "I"))))
  l1 <- DAISIE_A_loglik_element(pars, el, 3)
  expect_true(is.finite(l1))
  expect_equal(DAISIE_A_loglik_element(pars, el, 3, methode = "lsoda"), l1,
               tolerance = 1e-6)
})
