test_that("the recolonisation setting defaults to FALSE and is restored", {
  withr::local_options(DAISIE.recol_dd = NULL)
  expect_false(DAISIE_recol_dd())
  f <- function(recol_dd) {
    DAISIE:::local_recol_dd(recol_dd)
    DAISIE_recol_dd()
  }
  expect_true(f(TRUE))
  expect_false(DAISIE_recol_dd())
  withr::local_options(DAISIE.recol_dd = TRUE)
  expect_true(DAISIE_recol_dd())
  expect_false(f(FALSE))
  expect_true(DAISIE_recol_dd())
  expect_error(f(NA))
})

test_that("colonisation rate: recolonisation is not diversity dependent", {
  withr::local_options(DAISIE.recol_dd = FALSE)
  # 3 mainland species, 1 present as non-endemic, 2 species on the island
  pc <- 0.5 * (1 - 2 / 4)
  expect_equal(DAISIE:::get_immig_rate(gam = 0.5, num_spec = 2, K = 4,
                                       mainland_n = 3, num_immigrants = 1),
               2 * pc + 1 * 0.5)
  # with K reached, only recolonisation remains
  expect_equal(DAISIE:::get_immig_rate(gam = 0.5, num_spec = 4, K = 4,
                                       mainland_n = 3, num_immigrants = 1),
               0.5)
  withr::local_options(DAISIE.recol_dd = TRUE)
  expect_equal(DAISIE:::get_immig_rate(gam = 0.5, num_spec = 2, K = 4,
                                       mainland_n = 3, num_immigrants = 1),
               3 * pc)
  expect_null(DAISIE:::get_recol_prob(0.5, 1, 1))
})

test_that("colonist choice follows the recolonisation probability", {
  withr::local_options(DAISIE.recol_dd = FALSE)
  island_spec <- matrix(c("2", "2", "1.5", "I", NA, NA, NA), nrow = 1)
  set.seed(1)
  col <- replicate(2000, DAISIE:::choose_colonist(1:4, island_spec, 0.7))
  expect_equal(mean(col == 2), 0.7, tolerance = 0.05)
  expect_true(all(col %in% 1:4))
  # NULL: all mainland species equally likely (previous behaviour)
  set.seed(2)
  col <- replicate(2000, DAISIE:::choose_colonist(1:4, island_spec, NULL))
  expect_equal(mean(col == 2), 0.25, tolerance = 0.05)
})

test_that("CS likelihood: recol_dd = TRUE reproduces DAISIE 5.0.3", {
  pars1 <- c(0.6, 0.3, 8, 0.6, 0.4)
  f <- function(b, st, miss = 0, methode = "lsodes") {
    DAISIE:::DAISIE_loglik_CS_choice(
      pars1 = pars1, pars2 = c(100, 11, 0, 0), brts = b, stac = st,
      missnumspec = miss, methode = methode,
      abstolint = 1e-16, reltolint = 1e-12)
  }
  withr::local_options(DAISIE.recol_dd = TRUE)
  expect_equal(f(c(4, 2.5), 4), -5.1485895816, tolerance = 1e-8)
  expect_equal(f(c(4, 3, 1.2), 2, 1), -4.8964057017, tolerance = 1e-8)
  expect_equal(f(c(4, 2), 1), -1.9430826406, tolerance = 1e-8)
  expect_equal(f(c(4, 3.9, 1), 8), -3.3427584633, tolerance = 1e-8)
  withr::local_options(DAISIE.recol_dd = FALSE)
  expect_equal(f(c(4, 2.5), 4), -5.2992323334, tolerance = 1e-8)
  expect_equal(f(c(4, 3, 1.2), 2, 1), -4.9170633981, tolerance = 1e-8)
  expect_equal(f(c(4, 2), 1), -1.9353615013, tolerance = 1e-8)
  expect_equal(f(c(4, 3.9, 1), 8), -3.4198090124, tolerance = 1e-8)
  # R, Fortran and C++ right-hand sides agree
  for (st in c(4, 1, 8)) {
    b <- switch(as.character(st), "4" = c(4, 2.5), "1" = c(4, 2),
                "8" = c(4, 3.9, 1))
    l1 <- f(b, st, methode = "lsodes")
    expect_equal(f(b, st, methode = "deSolve_R::lsodes"), l1, tolerance = 1e-7)
    expect_equal(f(b, st, methode = "odeint::runge_kutta_fehlberg78"), l1,
                 tolerance = 1e-7)
  }
})

test_that("CS likelihood: the setting does not matter without diversity dependence of colonisation", {
  pars1 <- c(0.6, 0.3, 8, 0.6, 0.4)
  f <- function(ddm, K = 8) {
    DAISIE:::DAISIE_loglik_CS_choice(
      pars1 = c(0.6, 0.3, K, 0.6, 0.4), pars2 = c(100, ddm, 0, 0),
      brts = c(4, 2.5), stac = 4, missnumspec = 0, methode = "lsodes")
  }
  for (ddm in c(0, 1)) {
    l_new <- withr::with_options(list(DAISIE.recol_dd = FALSE), f(ddm))
    l_old <- withr::with_options(list(DAISIE.recol_dd = TRUE), f(ddm))
    expect_equal(l_new, l_old)
  }
  l_new <- withr::with_options(list(DAISIE.recol_dd = FALSE), f(11, Inf))
  l_old <- withr::with_options(list(DAISIE.recol_dd = TRUE), f(11, Inf))
  expect_equal(l_new, l_old)
})

test_that("CS likelihood: a maximum colonisation time integrates the known-time likelihood", {
  withr::local_options(DAISIE.recol_dd = FALSE)
  pars1 <- c(0.6, 0.5, 3, 1.5, 0.5)
  L <- function(b, stac) {
    exp(DAISIE:::DAISIE_loglik_CS_choice(
      pars1 = pars1, pars2 = c(100, 11, 0, 0), brts = b, stac = stac,
      missnumspec = 0, methode = "odeint::runge_kutta_fehlberg78",
      abstolint = 1e-14, reltolint = 1e-12))
  }
  # Gauss-Legendre quadrature on (0, 1.3)
  nodes <- c(-0.9739065285, -0.8650633667, -0.6794095683, -0.4333953941,
             -0.1488743390, 0.1488743390, 0.4333953941, 0.6794095683,
             0.8650633667, 0.9739065285)
  weights <- c(0.0666713443, 0.1494513492, 0.2190863625, 0.2692667193,
               0.2955242247, 0.2955242247, 0.2692667193, 0.2190863625,
               0.1494513492, 0.0666713443)
  tmax <- 1.3
  x <- tmax / 2 * nodes + tmax / 2
  integral <- tmax / 2 * sum(weights * vapply(x, function(t) L(c(2, t), 4), 0))
  expect_equal(L(c(2, tmax), 1), integral, tolerance = 1e-6)
})

test_that("IW likelihood: recol_dd = TRUE reproduces DAISIE 5.0.3", {
  dl <- list(list(island_age = 2, not_present = 1),
             list(colonist_name = "x", branching_times = c(2, 0.7), stac = 4,
                  missing_species = 0))
  f <- function(recol_dd, K = 3) {
    DAISIE_loglik_IW(pars1 = c(0.6, 0.5, K, 1.5, 0.5), pars2 = c(30, 11, 0, 0),
                     datalist = dl, methode = "ode45", recol_dd = recol_dd)
  }
  expect_equal(f(TRUE), -3.0641712278, tolerance = 1e-6)
  expect_equal(f(FALSE), -3.1754568395, tolerance = 1e-6)
  expect_equal(f(TRUE, Inf), f(FALSE, Inf))
})

test_that("simulation: recolonisation resets colonisation times at a carrying capacity of one", {
  # With K = 1 a non-endemic blocks all colonisation if recolonisation is
  # diversity dependent, so its colonisation time is its first arrival; if
  # not, it is recolonised at rate gam and its colonisation time is recent.
  sim_col_age <- function(recol_dd, use_rcpp) {
    set.seed(4)
    sims <- DAISIE_sim_cr(time = 3, M = 1, pars = c(0, 0.05, 1, 2, 0),
                          replicates = 200, divdepmodel = "CS",
                          plot_sims = FALSE, verbose = FALSE,
                          use_rcpp = use_rcpp, recol_dd = recol_dd)
    ages <- vapply(sims, function(r) {
      if (length(r) < 2 || r[[2]]$stac != 4) return(NA_real_)
      r[[2]]$branching_times[2]
    }, 0)
    mean(ages, na.rm = TRUE)
  }
  for (use_rcpp in c(FALSE, TRUE)) {
    expect_gt(sim_col_age(TRUE, use_rcpp), 2)
    expect_lt(sim_col_age(FALSE, use_rcpp), 1)
  }
  expect_false(DAISIE_recol_dd())
})
