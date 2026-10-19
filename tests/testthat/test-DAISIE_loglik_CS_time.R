# With constant area (no island ontogeny, no sea level) the time-dependent
# likelihood must be identical to the time-independent one with linear
# diversity-dependence in cladogenesis and immigration (ddep = 11).
test_that("time-dependent likelihood equals time-independent one for constant area", {
  island_age <- 10
  lx <- 30
  pars_static <- c(0.5, 0.3, 20, 0.05, 0.8) # lac, mu, K, gam, laa
  pars_time <- c(
    pars_static,
    0, 0, # d, x
    1, 1, 0, island_age, 0, 0, 0, # area pars: max, current, peak, age, ...
    0, 0, island_age, 1 # ontogeny, sea level, total time, peak
  )
  cases <- list(
    list(stac = 1, brts = c(island_age, island_age, 0)),
    list(stac = 5, brts = c(island_age, island_age, 0)),
    list(stac = 6, brts = c(island_age, island_age, 3, 0)),
    list(stac = 7, brts = c(island_age, island_age, 3, 0)),
    list(stac = 8, brts = c(island_age, island_age, 4, 0)),
    list(stac = 9, brts = c(island_age, island_age, 4, 0))
  )
  for (case in cases) {
    loglik_static <- DAISIE_loglik_CS_M1(
      pars1 = pars_static,
      pars2 = c(lx, 11, 0, 0, NA, NA),
      brts = case$brts,
      stac = case$stac,
      missnumspec = 0,
      methode = "lsodes"
    )
    loglik_time <- DAISIE_loglik_CS_M1(
      pars1 = pars_time,
      pars2 = c(lx, 11, 0, 0, 0, 0),
      brts = case$brts,
      stac = case$stac,
      missnumspec = 0,
      methode = "lsodes"
    )
    testthat::expect_equal(
      loglik_time, loglik_static, tolerance = 1e-6,
      info = paste("stac =", case$stac)
    )
  }
})

# Helper: time-dependent parameter vector with constant area
cs_time_const_area_pars <- function(pars_static, island_age) {
  c(pars_static, 0, 0, 1, 1, 0, island_age, 0, 0, 0, 0, 0, island_age, 1)
}

# The three time-dependent right-hand sides must reproduce the static ones
# derivative by derivative when the area is constant.
test_that("time-dependent rhs functions equal static rhs functions at constant area", {
  testthat::skip_on_cran()
  island_age <- 10
  lx <- 7
  pars_static <- c(0.5, 0.3, 20, 0.05, 0.8)
  pars_time <- cs_time_const_area_pars(pars_static, island_age)
  funs <- list(
    list(stat = DAISIE_loglik_rhs,  time = DAISIE_loglik_rhs_time,  len = 2 * lx + 1),
    list(stat = DAISIE_loglik_rhs1, time = DAISIE_loglik_rhs_time1, len = 4 * lx),
    list(stat = DAISIE_loglik_rhs2, time = DAISIE_loglik_rhs_time2, len = 3 * lx)
  )
  set.seed(42)
  for (kk in 0:1) {
    parsvec_static <- DAISIE_loglik_rhs_precomp(c(pars_static, kk, 11), lx)
    parsvec_time <- c(pars_time, kk, 11)
    for (f in funs) {
      for (rep in 1:3) {
        x <- stats::runif(f$len)
        for (tt in c(-9, -4, -0.5)) {
          d_static <- f$stat(0, x, parsvec_static)[[1]]
          d_time <- f$time(tt, x, parsvec_time)[[1]]
          testthat::expect_equal(
            d_time, d_static, tolerance = 1e-10,
            info = paste("kk =", kk, "len =", f$len, "t =", tt)
          )
        }
      }
    }
  }
})

# Constant-area equivalence for every stac, with missing species and for
# colonisation times younger than the island.
test_that("time-dependent likelihood equals static one for all stac values", {
  testthat::skip_on_cran()
  island_age <- 10
  pars_static <- c(0.5, 0.3, 20, 0.05, 0.8)
  pars_time <- cs_time_const_area_pars(pars_static, island_age)
  cases <- list(
    list(stac = 0, brts = c(island_age, 0), miss = 0),
    list(stac = 2, brts = c(island_age, 4, 0), miss = 0),
    list(stac = 3, brts = c(island_age, 6, 4, 0), miss = 0),
    list(stac = 4, brts = c(island_age, 6, 4, 0), miss = 0),
    list(stac = 1, brts = c(island_age, 6, 0), miss = 0),
    list(stac = 5, brts = c(island_age, 6, 0), miss = 0),
    list(stac = 6, brts = c(island_age, 6, 3, 0), miss = 0),
    list(stac = 7, brts = c(island_age, 6, 3, 0), miss = 0),
    list(stac = 8, brts = c(island_age, 6, 3, 0), miss = 0),
    list(stac = 9, brts = c(island_age, 6, 3, 0), miss = 0),
    list(stac = 2, brts = c(island_age, 4, 0), miss = 2),
    list(stac = 1, brts = c(island_age, island_age, 0), miss = 3),
    list(stac = 5, brts = c(island_age, island_age, 0), miss = 3),
    list(stac = 6, brts = c(island_age, island_age, 3, 0), miss = 2),
    list(stac = 9, brts = c(island_age, 6, 3, 0), miss = 2)
  )
  for (lx in c(25, 40)) {
    for (case in cases) {
      loglik_static <- DAISIE_loglik_CS_M1(
        pars1 = pars_static, pars2 = c(lx, 11, 0, 0, NA, NA),
        brts = case$brts, stac = case$stac, missnumspec = case$miss,
        methode = "lsodes"
      )
      loglik_time <- DAISIE_loglik_CS_M1(
        pars1 = pars_time, pars2 = c(lx, 11, 0, 0, 0, 0),
        brts = case$brts, stac = case$stac, missnumspec = case$miss,
        methode = "lsodes"
      )
      testthat::expect_equal(
        loglik_time, loglik_static, tolerance = 1e-6,
        info = paste("stac =", case$stac, "miss =", case$miss, "lx =", lx,
                     "brts =", paste(case$brts, collapse = ","))
      )
    }
  }
})

# Parameters with beta-function island ontogeny (area changes through time).
# Layout: lac, mu, K, gam, laa, d, x, area_pars (max, current, peak, age,
# sea-level amplitude, frequency, angle), ontogeny, sea level, total time, peak
beta_ontogeny_pars <- function(island_age) {
  c(0.5, 0.3, 20, 0.05, 0.8, 0.1, 0.1,
    1, 0.5, 0.5, island_age, 0, 0, 0,
    1, 0, island_age, 0.5)
}

# The ODE solvers must agree on the time-dependent likelihood.
test_that("time-dependent likelihood is solver independent", {
  testthat::skip_on_cran()
  island_age <- 10
  pars_time <- beta_ontogeny_pars(island_age)
  cases <- list(
    list(stac = 2, brts = c(island_age, 4, 0)),
    list(stac = 1, brts = c(island_age, island_age, 0)),
    list(stac = 6, brts = c(island_age, island_age, 3, 0)),
    list(stac = 3, brts = c(island_age, 6, 4, 0))
  )
  ll <- function(case, methode) {
    DAISIE_loglik_CS_M1(
      pars1 = pars_time, pars2 = c(30, 11, 0, 0, 1, 0),
      brts = case$brts, stac = case$stac, missnumspec = 0,
      methode = methode, reltolint = 1e-10, abstolint = 1e-14
    )
  }
  for (case in cases) {
    ref <- ll(case, "lsodes")
    testthat::expect_true(is.finite(ref))
    testthat::expect_equal(ll(case, "lsoda"), ref, tolerance = 1e-6,
                           info = paste("lsoda, stac", case$stac))
    testthat::expect_equal(ll(case, "ode45"), ref, tolerance = 1e-6,
                           info = paste("ode45, stac", case$stac))
  }
})

# DAISIE_integrate must pick the time-dependent path only for long pars.
test_that("DAISIE_integrate dispatches on the length of pars", {
  testthat::skip_on_cran()
  island_age <- 10
  lx <- 12
  pars_static <- c(0.5, 0.3, 20, 0.05, 0.8)
  init <- c(1, rep(0, 2 * lx))
  y_const <- DAISIE:::DAISIE_integrate(
    init, c(-island_age, 0), DAISIE_loglik_rhs, c(pars_static, 0, 11),
    rtol = 1e-10, atol = 1e-16, method = "lsodes"
  )
  y_time <- DAISIE:::DAISIE_integrate(
    init, c(-island_age, 0), DAISIE_loglik_rhs, # callers pass the static rhs
    c(cs_time_const_area_pars(pars_static, island_age), 0, 11),
    rtol = 1e-10, atol = 1e-16, method = "lsodes"
  )
  testthat::expect_equal(length(y_const), 2 * lx + 1)
  testthat::expect_equal(y_time, y_const, tolerance = 1e-6)
  # the dispatch in DAISIE_integrate(_time) greps the second statement of the
  # static rhs functions for these literals
  for (i in 0:2) {
    f <- list(DAISIE_loglik_rhs, DAISIE_loglik_rhs1, DAISIE_loglik_rhs2)[[i + 1]]
    testthat::expect_match(as.character(body(f)[2]), paste0("rhs <- ", i))
  }
})

# With a truly changing area, the time-dependent solution must be the limit of
# many short constant-rate steps in which the area is frozen at the midpoint.
test_that("time-dependent ODE is the limit of piecewise constant-area steps", {
  testthat::skip_on_cran()
  island_age <- 10
  lx <- 15
  pars_time <- beta_ontogeny_pars(island_age)
  lac <- pars_time[1]; mu <- pars_time[2]; K <- pars_time[3]
  gam <- pars_time[4]; laa <- pars_time[5]; d <- pars_time[6]; xh <- pars_time[7]
  init <- c(1, rep(0, 2 * lx))
  y_time <- DAISIE:::DAISIE_integrate(
    init, c(-island_age, 0), DAISIE_loglik_rhs, c(pars_time, 0, 11),
    rtol = 1e-12, atol = 1e-16, method = "lsodes"
  )
  stepwise <- function(nsteps) {
    edges <- seq(-island_age, 0, length.out = nsteps + 1)
    y <- init
    for (i in seq_len(nsteps)) {
      mid <- mean(edges[i:(i + 1)])
      A <- DAISIE:::island_area_vector(
        timeval = pars_time[17] - abs(mid), area_pars = pars_time[8:14],
        island_ontogeny = pars_time[15], sea_level = pars_time[16],
        total_time = pars_time[17], peak = pars_time[18]
      )
      # linear diversity dependence: per-capita rates scale with the area
      pars_step <- c(lac * A^d, mu * A^(-xh), K * A, gam, laa, 0, 11)
      y <- DAISIE:::DAISIE_integrate(
        y, edges[i:(i + 1)], DAISIE_loglik_rhs, pars_step,
        rtol = 1e-12, atol = 1e-16, method = "lsodes"
      )
    }
    y
  }
  err <- vapply(c(5, 20, 80), function(n) max(abs(stepwise(n) - y_time)),
                numeric(1))
  testthat::expect_true(all(diff(err) < 0), info = paste(err, collapse = " "))
  testthat::expect_lt(err[3], 1e-3)
  # the area really changes, so this is not the trivial constant-area case
  testthat::expect_gt(max(abs(y_time - DAISIE:::DAISIE_integrate(
    init, c(-island_age, 0), DAISIE_loglik_rhs,
    c(cs_time_const_area_pars(pars_time[1:5], island_age), 0, 11),
    rtol = 1e-12, atol = 1e-16, method = "lsodes"))), 1e-3)
})

# Regression values for a beta-ontogeny island (computed with the corrected
# state-vector layout and time direction; they guard against silent changes in the ODE system).
test_that("time-dependent likelihood reproduces stored reference values", {
  testthat::skip_on_cran()
  age <- 5
  total <- 10
  peak <- calc_peak(age, create_area_pars(1000, 500, 0.1, total, 0, 0, 0))
  pars <- c(0.5, 0.3, 20, 0.15, 0.8, 0.1, 0.1, 1000, 500, 0.1, total, 0, 0, 0,
            1, 0, age, peak)
  ll <- function(brts, stac) {
    DAISIE_loglik_CS_M1(
      pars1 = pars, pars2 = c(40, 11, 0, 0, 1, 0), brts = brts, stac = stac,
      missnumspec = 0, methode = "lsodes"
    )
  }
  testthat::expect_equal(ll(c(age, 3, 0), 2), -5.2936849440, tolerance = 1e-6)
  testthat::expect_equal(ll(c(age, age, 0), 1), -3.1206793922, tolerance = 1e-6)
  testthat::expect_equal(ll(c(age, age, 2, 0), 6), -6.2143412781, tolerance = 1e-6)
})

# The probability of an empty island predicted by the likelihood (stac 0) must
# match the fraction of empty islands in simulations. The area-dependence of
# extinction and cladogenesis is strong here and the area changes a lot, so
# this test detects a wrong direction of time in the area function: area is
# zero at the emergence of the island and equals current_area at the present.
test_that("simulated fraction of empty islands matches the likelihood", {
  testthat::skip_on_cran()
  age <- 5
  total <- 10
  replicates <- 300
  pars <- c(0.5, 0.3, 5, 0.4, 0.5)
  area_pars <- create_area_pars(2, 1, 0.3, total, 0, 0, 0)
  set.seed(3)
  sim <- DAISIE_sim_time_dep(
    time = age, M = 1, pars = pars, replicates = replicates,
    island_ontogeny = "beta", sea_level = "const",
    area_pars = area_pars,
    hyper_pars = create_hyper_pars(d = 0.3, x = 1.2),
    plot_sims = FALSE, verbose = FALSE
  )
  frac_empty <- mean(vapply(sim, function(r) r[[1]]$not_present == 1, logical(1)))
  peak <- calc_peak(age, area_pars)
  pars_time <- c(pars, 0.3, 1.2, 2, 1, 0.3, total, 0, 0, 0, 1, 0, age, peak)
  p_empty <- exp(DAISIE_loglik_CS_M1(
    pars1 = pars_time, pars2 = c(40, 11, 0, 0, 1, 0), brts = c(age, 0),
    stac = 0, missnumspec = 0, methode = "odeint::runge_kutta_dopri5"
  ))
  se <- sqrt(p_empty * (1 - p_empty) / replicates)
  testthat::expect_lt(abs(frac_empty - p_empty), 4 * se)
})
