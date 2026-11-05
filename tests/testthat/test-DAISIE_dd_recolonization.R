test_that("dd_recolonization is an element of CS_version, FALSE by default", {
  expect_false(create_CS_version()$dd_recolonization)
  expect_true(create_CS_version(dd_recolonization = TRUE)$dd_recolonization)
  expect_true(create_CS_version(model = 2, relaxed_par = "cladogenesis",
                                par_sd = 1,
                                dd_recolonization = TRUE)$dd_recolonization)
  expect_error(create_CS_version(dd_recolonization = NA), "dd_recolonization")
  expect_error(create_CS_version(dd_recolonization = "yes"),
               "dd_recolonization")
  expect_false(get_dd_recolonization(1))
  expect_false(get_dd_recolonization(list(model = 1)))
  expect_true(get_dd_recolonization(list(model = 1, dd_recolonization = TRUE)))
  # no separate argument (and no option) any more
  expect_false("recol_dd" %in% names(formals(DAISIE_loglik_CS)))
  expect_false("dd_recolonization" %in% names(formals(DAISIE_loglik_CS)))
  expect_false("dd_recolonization" %in% names(formals(DAISIE_ML_CS)))
  expect_false(exists("DAISIE_recol_dd", envir = asNamespace("DAISIE")))
  expect_false(formals(DAISIE_loglik_IW)$dd_recolonization)
  expect_false(formals(DAISIE_sim_cr)$dd_recolonization)
  expect_error(DAISIE_sim_cr(time = 1, M = 1, pars = c(1, 1, 1, 1, 1),
                             replicates = 1, plot_sims = FALSE,
                             verbose = FALSE, dd_recolonization = NA),
               "dd_recolonization")
})

test_that("colonization rate: recolonization is not diversity dependent", {
  # 3 mainland species, 1 present as non-endemic, 2 species on the island
  pc <- 0.5 * (1 - 2 / 4)
  expect_equal(get_immig_rate(gam = 0.5, num_spec = 2, K = 4,
                              mainland_n = 3, num_immigrants = 1),
               2 * pc + 1 * 0.5)
  # with K reached, only recolonization remains
  expect_equal(get_immig_rate(gam = 0.5, num_spec = 4, K = 4,
                              mainland_n = 3, num_immigrants = 1),
               0.5)
  expect_equal(get_immig_rate(gam = 0.5, num_spec = 2, K = 4,
                              mainland_n = 3, num_immigrants = 1,
                              dd_recolonization = TRUE),
               3 * pc)
  expect_null(get_recol_prob(0.5, 1, 1, dd_recolonization = TRUE))
  expect_equal(get_recol_prob(0.5, 1, 1), 0.5)
})

test_that("colonist choice follows the recolonization probability", {
  island_spec <- matrix(c("2", "2", "1.5", "I", NA, NA, NA), nrow = 1)
  set.seed(1)
  col <- replicate(2000, choose_colonist(1:4, island_spec, 0.7))
  expect_equal(mean(col == 2), 0.7, tolerance = 0.05)
  expect_true(all(col %in% 1:4))
  # NULL: all mainland species equally likely (as with dd_recolonization = TRUE)
  set.seed(2)
  col <- replicate(2000, choose_colonist(1:4, island_spec, NULL))
  expect_equal(mean(col == 2), 0.25, tolerance = 0.05)
})

test_that("CS likelihood: dd_recolonization = TRUE reproduces DAISIE 5.0.3", {
  pars1 <- c(0.6, 0.3, 8, 0.6, 0.4)
  f <- function(b, st, miss = 0, methode = "lsodes", dd_recolonization) {
    DAISIE_loglik_CS_choice(
      pars1 = pars1, pars2 = c(100, 11, 0, 0), brts = b, stac = st,
      missnumspec = miss, methode = methode,
      CS_version = create_CS_version(dd_recolonization = dd_recolonization),
      abstolint = 1e-16, reltolint = 1e-12)
  }
  # values of DAISIE 5.0.3 (develop 6e3c1ab)
  expect_equal(f(c(4, 2.5), 4, dd_recolonization = TRUE), -5.1485895816,
               tolerance = 1e-8)
  expect_equal(f(c(4, 3, 1.2), 2, 1, dd_recolonization = TRUE), -4.8964057017,
               tolerance = 1e-8)
  expect_equal(f(c(4, 2), 1, dd_recolonization = TRUE), -1.9430826406,
               tolerance = 1e-8)
  expect_equal(f(c(4, 3.9, 1), 8, dd_recolonization = TRUE), -3.3427584633,
               tolerance = 1e-8)
  expect_equal(f(c(4, 2.5), 4, dd_recolonization = FALSE), -5.2992323334,
               tolerance = 1e-8)
  expect_equal(f(c(4, 3, 1.2), 2, 1, dd_recolonization = FALSE), -4.9170633981,
               tolerance = 1e-8)
  expect_equal(f(c(4, 2), 1, dd_recolonization = FALSE), -1.9353615013,
               tolerance = 1e-8)
  expect_equal(f(c(4, 3.9, 1), 8, dd_recolonization = FALSE), -3.4198090124,
               tolerance = 1e-8)
  # R, Fortran and C++ right-hand sides agree
  for (dd_recolonization in c(FALSE, TRUE)) {
    for (st in c(4, 1, 8)) {
      b <- switch(as.character(st), "4" = c(4, 2.5), "1" = c(4, 2),
                  "8" = c(4, 3.9, 1))
      l1 <- f(b, st, methode = "lsodes", dd_recolonization = dd_recolonization)
      expect_equal(f(b, st, methode = "deSolve_R::lsodes",
                     dd_recolonization = dd_recolonization), l1,
                   tolerance = 1e-7)
      expect_equal(f(b, st, methode = "odeint::runge_kutta_fehlberg78",
                     dd_recolonization = dd_recolonization), l1,
                   tolerance = 1e-7)
    }
  }
  # the default is FALSE, also when CS_version is a number
  expect_equal(DAISIE_loglik_CS_choice(
    pars1 = pars1, pars2 = c(100, 11, 0, 0), brts = c(4, 2.5), stac = 4,
    missnumspec = 0, methode = "lsodes", CS_version = 1,
    abstolint = 1e-16, reltolint = 1e-12),
    f(c(4, 2.5), 4, dd_recolonization = FALSE))
})

test_that("CS likelihood: the setting reaches all clades and the ML", {
  utils::data(Galapagos_datalist, package = "DAISIE")
  pars1 <- c(2.5, 2.7, 20, 0.009, 1.01)
  pars2 <- c(40, 11, 0, 0)
  l_old <- DAISIE_loglik_CS(pars1, pars2, Galapagos_datalist,
                            methode = "odeint::runge_kutta_fehlberg78",
                            CS_version = create_CS_version(
                              dd_recolonization = TRUE))
  l_new <- DAISIE_loglik_CS(pars1, pars2, Galapagos_datalist,
                            methode = "odeint::runge_kutta_fehlberg78")
  l_new2 <- DAISIE_loglik_CS(pars1, pars2, Galapagos_datalist,
                             methode = "odeint::runge_kutta_fehlberg78",
                             CS_version = create_CS_version(
                               dd_recolonization = FALSE))
  expect_identical(l_new, l_new2)
  # sum over the clades
  per_clade <- function(dd_recolonization) {
    logp0 <- DAISIE_loglik_CS_choice(
      pars1 = pars1, pars2 = pars2, brts = Galapagos_datalist[[1]]$island_age,
      stac = 0, missnumspec = 0, methode = "odeint::runge_kutta_fehlberg78",
      CS_version = create_CS_version(dd_recolonization = dd_recolonization))
    l <- Galapagos_datalist[[1]]$not_present * logp0
    for (i in 2:length(Galapagos_datalist)) {
      l <- l + DAISIE_loglik_CS_choice(
        pars1 = pars1, pars2 = pars2,
        brts = Galapagos_datalist[[i]]$branching_times,
        stac = Galapagos_datalist[[i]]$stac,
        missnumspec = Galapagos_datalist[[i]]$missing_species,
        methode = "odeint::runge_kutta_fehlberg78",
        CS_version = create_CS_version(dd_recolonization = dd_recolonization))
    }
    l
  }
  expect_equal(l_old, per_clade(TRUE), tolerance = 1e-10)
  expect_equal(l_new, per_clade(FALSE), tolerance = 1e-10)
  expect_gt(abs(l_new - l_old), 1e-6)
  # DAISIE_ML_CS passes CS_version on to the likelihood: the maximum
  # loglikelihood is the loglikelihood (with the same setting) at the estimate
  skip_on_cran()
  for (dd_recolonization in c(FALSE, TRUE)) {
    cs_version <- create_CS_version(dd_recolonization = dd_recolonization)
    ml <- DAISIE_ML_CS(
      datalist = Galapagos_datalist, initparsopt = pars1[4], idparsopt = 4,
      parsfix = pars1[-4], idparsfix = c(1, 2, 3, 5), ddmodel = 11, res = 40,
      methode = "odeint::runge_kutta_fehlberg78", verbose = 0,
      CS_version = cs_version)
    expect_equal(ml$loglik,
                 DAISIE_loglik_CS(c(pars1[1:3], ml$gamma, pars1[5]), pars2,
                                  Galapagos_datalist,
                                  methode = "odeint::runge_kutta_fehlberg78",
                                  CS_version = cs_version),
                 tolerance = 1e-8)
    other <- DAISIE_loglik_CS(c(pars1[1:3], ml$gamma, pars1[5]), pars2,
                              Galapagos_datalist,
                              methode = "odeint::runge_kutta_fehlberg78",
                              CS_version = create_CS_version(
                                dd_recolonization = !dd_recolonization))
    expect_gt(abs(ml$loglik - other), 1e-6)
  }
})

test_that("CS likelihood: the setting does not matter without diversity dependence of colonization", {
  f <- function(ddm, K = 8, dd_recolonization) {
    DAISIE_loglik_CS_choice(
      pars1 = c(0.6, 0.3, K, 0.6, 0.4), pars2 = c(100, ddm, 0, 0),
      brts = c(4, 2.5), stac = 4, missnumspec = 0, methode = "lsodes",
      CS_version = create_CS_version(dd_recolonization = dd_recolonization))
  }
  for (ddm in c(0, 1)) {
    expect_equal(f(ddm, dd_recolonization = FALSE),
                 f(ddm, dd_recolonization = TRUE))
  }
  expect_equal(f(11, Inf, dd_recolonization = FALSE),
               f(11, Inf, dd_recolonization = TRUE))
})

test_that("CS likelihood: a maximum colonization time integrates the known-time likelihood", {
  pars1 <- c(0.6, 0.5, 3, 1.5, 0.5)
  for (dd_recolonization in c(FALSE, TRUE)) {
    L <- function(b, stac) {
      exp(DAISIE_loglik_CS_choice(
        pars1 = pars1, pars2 = c(100, 11, 0, 0), brts = b, stac = stac,
        missnumspec = 0, methode = "odeint::runge_kutta_fehlberg78",
        CS_version = create_CS_version(dd_recolonization = dd_recolonization),
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
  }
})

test_that("IW likelihood: dd_recolonization = TRUE reproduces DAISIE 5.0.3", {
  dl <- list(list(island_age = 2, not_present = 1),
             list(colonist_name = "x", branching_times = c(2, 0.7), stac = 4,
                  missing_species = 0))
  f <- function(dd_recolonization, K = 3) {
    DAISIE_loglik_IW(pars1 = c(0.6, 0.5, K, 1.5, 0.5), pars2 = c(30, 11, 0, 0),
                     datalist = dl, methode = "ode45",
                     dd_recolonization = dd_recolonization)
  }
  expect_equal(f(TRUE), -3.0641712278, tolerance = 1e-6)
  expect_equal(f(FALSE), -3.1754568395, tolerance = 1e-6)
  expect_equal(f(TRUE, Inf), f(FALSE, Inf))
  # the default is FALSE
  expect_identical(
    DAISIE_loglik_IW(pars1 = c(0.6, 0.5, 3, 1.5, 0.5), pars2 = c(30, 11, 0, 0),
                     datalist = dl, methode = "ode45"),
    f(FALSE))
})

test_that("simulation: recolonization resets colonization times at a carrying capacity of one", {
  # With K = 1 a non-endemic blocks all colonization if recolonization is
  # diversity dependent, so its colonization time is its first arrival; if
  # not, it is recolonized at rate gam and its colonization time is recent.
  sim_col_age <- function(dd_recolonization, use_rcpp) {
    set.seed(4)
    sims <- DAISIE_sim_cr(time = 3, M = 1, pars = c(0, 0.05, 1, 2, 0),
                          replicates = 200, divdepmodel = "CS",
                          plot_sims = FALSE, verbose = FALSE,
                          use_rcpp = use_rcpp,
                          dd_recolonization = dd_recolonization)
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
})
