# Does the time-dependent (island area) model reduce to the time-constant model
# when the exponents of the area dependence, d (cladogenesis) and x
# (extinction), go to zero?
#
# In the time-dependent model the per-capita rates are
#   cladogenesis  lac * A^d * (1 - n / (K * A))
#   extinction    mu  * A^(-x)
#   immigration   gam * (1 - n / (K * A))
# so for d = x = 0 the area still enters through the carrying capacity K * A(t).
# The time-constant model is therefore only recovered at d = x = 0 when K = Inf
# (no diversity dependence) or when the area is constant. All tests below use
# K = Inf and an area that changes strongly through time (beta-function island
# ontogeny, with and without sea level).

age <- 5          # island age
lx <- 30

# area parameters: max area, current area, proportional peak time, total island
# age, sea level amplitude, sea level frequency, island gradient angle
area_scenarios <- list(
  beta     = list(ap = create_area_pars(1000, 500, 0.3, 10, 0, 0, 0),  sea = 0),
  beta_sea = list(ap = create_area_pars(1000, 500, 0.1, 10, 1, 2, 10), sea = 1),
  # strongly peaked beta function (shape parameter ~114): the area underflows to
  # exactly 0 during the first ~1e-4 time units after island emergence
  beta_low = list(ap = create_area_pars(1000, 100, 0.6, 10, 0, 0, 0),  sea = 0)
)
ap_vec <- function(ap) {
  unlist(ap)[c("max_area", "current_area", "proportional_peak_t",
               "total_island_age", "sea_level_amplitude",
               "sea_level_frequency", "island_gradient_angle")]
}

# loglik of one clade, time-constant and time-dependent version
ll_clade_const <- function(p5, brts, stac, miss, methode) {
  DAISIE_loglik_CS_M1(
    pars1 = p5, pars2 = c(lx, 11, 0, 0, NA, NA), brts = brts, stac = stac,
    missnumspec = miss, methode = methode
  )
}
ll_clade_time <- function(p5, d, x, sc, brts, stac, miss, methode) {
  peak <- calc_peak(age, sc$ap)
  DAISIE_loglik_CS_M1(
    pars1 = c(p5, d, x, ap_vec(sc$ap), 1, sc$sea, age, peak),
    pars2 = c(lx, 11, 0, 0, 1, sc$sea), brts = brts, stac = stac,
    missnumspec = miss, methode = methode
  )
}
# loglik of a whole data set
ll_all_const <- function(p5, dl, methode) {
  DAISIE_loglik_all(pars1 = p5, pars2 = c(40, 11, 0, 0), datalist = dl,
                    methode = methode)
}
ll_all_time <- function(p5, d, x, sc, dl, methode) {
  DAISIE_loglik_all(
    pars1 = c(p5, d, x, ap_vec(sc$ap)), pars2 = c(40, 11, 0, 0, 1, sc$sea),
    datalist = dl, methode = methode
  )
}

clade_cases <- list(
  list(stac = 0, brts = c(age, 0),         miss = 0),
  list(stac = 2, brts = c(age, 3, 0),      miss = 0),
  list(stac = 2, brts = c(age, 3, 0),      miss = 2),
  list(stac = 3, brts = c(age, 4, 2, 0),   miss = 2),
  list(stac = 4, brts = c(age, 3, 0),      miss = 0),
  list(stac = 1, brts = c(age, age, 0),    miss = 3),
  list(stac = 5, brts = c(age, age, 0),    miss = 3),
  list(stac = 1, brts = c(age, 3.5, 0),    miss = 0), # max colonisation time
  list(stac = 5, brts = c(age, 3.5, 0),    miss = 2),
  list(stac = 6, brts = c(age, 4, 2, 0),   miss = 2),
  list(stac = 7, brts = c(age, age, 3, 0), miss = 0),
  list(stac = 8, brts = c(age, 4, 2, 0),   miss = 0),
  list(stac = 9, brts = c(age, 4, 2, 0),   miss = 2)
)
par_sets <- list(
  c(0.5, 0.3, Inf, 0.10, 0.4),
  c(1.5, 0.8, Inf, 0.05, 1.0)
)
solver <- "odeint::runge_kutta_dopri5"

test_that("time-dependent clade loglik equals time-constant one at d = x = 0, K = Inf", {
  for (sn in c("beta", "beta_sea")) {
    for (p5 in par_sets) {
      for (cs in clade_cases) {
        expect_equal(
          ll_clade_time(p5, 0, 0, area_scenarios[[sn]], cs$brts, cs$stac,
                        cs$miss, solver),
          ll_clade_const(p5, cs$brts, cs$stac, cs$miss, solver),
          tolerance = 1e-10,
          info = paste(sn, "stac =", cs$stac, "miss =", cs$miss)
        )
      }
    }
  }
})

test_that("time-dependent data set loglik equals time-constant one at d = x = 0, K = Inf", {
  skip_on_cran()
  p5 <- c(0.5, 0.3, Inf, 0.1, 0.4)
  set.seed(1)
  dl <- DAISIE_sim_cr(time = age, M = 40, pars = p5, replicates = 1,
                      plot_sims = FALSE, verbose = FALSE)[[1]]
  for (sn in c("beta", "beta_sea")) {
    expect_equal(
      ll_all_time(p5, 0, 0, area_scenarios[[sn]], dl, solver),
      ll_all_const(p5, dl, solver),
      tolerance = 1e-10, info = sn
    )
  }
})

test_that("time-dependent loglik is continuous in d and x around 0 (K = Inf)", {
  skip_on_cran()
  p5 <- c(0.5, 0.3, Inf, 0.1, 0.4)
  set.seed(1)
  dl <- DAISIE_sim_cr(time = age, M = 40, pars = p5, replicates = 1,
                      plot_sims = FALSE, verbose = FALSE)[[1]]
  for (sn in c("beta", "beta_sea")) {
    ll0 <- ll_all_const(p5, dl, solver)
    for (eps in c(1e-3, 1e-5, 1e-7, 1e-9)) {
      for (dx in list(c(eps, eps), c(eps, 0), c(0, eps))) {
        # the difference vanishes linearly in eps (here with a slope below 20)
        expect_lt(
          abs(ll_all_time(p5, dx[1], dx[2], area_scenarios[[sn]], dl, solver) - ll0),
          20 * eps + 1e-9
        )
      }
    }
  }
})

# The next two tests use an area that underflows to exactly 0 early on in the
# island's history. Then
# (1) K * A = Inf * 0 = NaN in the diversity dependence for K = Inf (and
#     n / (K * A) = 0 / 0 for n = 0 and finite K), which
#     pmax(0, NaN, na.rm = TRUE) turns into a cladogenesis and immigration rate
#     of 0 instead of lac * A^d and gam. The first test below fails because of
#     this (loglik off by ~1e-5 per clade, ~4e-4 for the data set).
# (2) A^(-x) = Inf for any x > 0, however small, which is cut off at the
#     extcutoff of 1e6, but equals mu for x = 0. This discontinuity in x is
#     currently masked by (1), because no species can be present when the
#     immigration rate is 0, so the second test passes; but it appears as soon
#     as only (1) is fixed (e.g. by guarding the NaN) and makes the likelihood
#     depend on the ODE solver.
# Both are fixed by flooring the area at a small positive number (e.g. 1e-300)
# in island_area_vector() (R) and in island_area_vector() of
# src/DAISIE_time_rates.h (C++).
test_that("equality at d = x = 0, K = Inf also holds when the area underflows to 0", {
  p5 <- c(0.5, 0.3, Inf, 0.1, 0.4)
  for (cs in clade_cases) {
    expect_equal(
      ll_clade_time(p5, 0, 0, area_scenarios$beta_low, cs$brts, cs$stac,
                    cs$miss, solver),
      ll_clade_const(p5, cs$brts, cs$stac, cs$miss, solver),
      tolerance = 1e-10,
      info = paste("stac =", cs$stac, "miss =", cs$miss)
    )
  }
})

test_that("loglik is continuous at x = 0 when the area underflows to 0", {
  skip_on_cran()
  p5 <- c(0.5, 0.3, Inf, 0.1, 0.4)
  set.seed(1)
  dl <- DAISIE_sim_cr(time = age, M = 40, pars = p5, replicates = 1,
                      plot_sims = FALSE, verbose = FALSE)[[1]]
  sc <- area_scenarios$beta_low
  ll0 <- ll_all_time(p5, 0, 0, sc, dl, "odeint::runge_kutta_cash_karp54")
  for (eps in c(1e-5, 1e-6, 1e-7, 1e-8)) {
    # the slope in x is ~40 here because |log(A)| is large
    expect_lt(
      abs(ll_all_time(p5, 0, eps, sc, dl, "odeint::runge_kutta_cash_karp54") - ll0),
      100 * eps + 1e-9
    )
  }
})

# Maximum likelihood ---------------------------------------------------------

test_that("ML of the time-dependent model at d = x = 0 and K = Inf equals ML of the time-constant model", {
  skip_on_cran()
  set.seed(1)
  dl <- DAISIE_sim_cr(time = age, M = 100, pars = c(0.4, 0.2, Inf, 0.04, 0.5),
                      replicates = 1, plot_sims = FALSE, verbose = FALSE)[[1]]
  start <- c(0.3, 0.3, 0.06, 0.3)
  tol <- c(1e-8, 1e-8, 1e-8)
  ml_const <- DAISIE_ML(
    datalist = dl, initparsopt = start, idparsopt = c(1, 2, 4, 5),
    parsfix = Inf, idparsfix = 3, ddmodel = 11, res = 40, cond = 0,
    tol = tol, methode = solver, optimmethod = "simplex", verbose = 0
  )
  ml_time <- function(sc, d, x) {
    DAISIE_ML(
      datalist = dl, initparsopt = start, idparsopt = c(1, 2, 4, 5),
      parsfix = c(Inf, d, x), idparsfix = c(3, 6, 7), island_ontogeny = 1,
      sea_level = sc$sea, area_pars = sc$ap, ddmodel = 11, res = 40, cond = 0,
      tol = tol, methode = solver, optimmethod = "simplex", verbose = 0
    )
  }
  pars <- c("lambda_c", "mu", "gamma", "lambda_a")
  for (sn in c("beta", "beta_sea")) {
    for (eps in c(0, 1e-6)) {
      fit <- ml_time(area_scenarios[[sn]], eps, eps)
      info <- paste(sn, "d = x =", eps)
      expect_equal(fit$conv, 0, info = info)
      expect_equal(fit$loglik, ml_const$loglik, tolerance = 1e-6, info = info)
      expect_equal(unlist(fit[pars]), unlist(ml_const[pars]), tolerance = 1e-3,
                   info = info)
    }
  }
})
