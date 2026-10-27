# Islands with sea level change but without island ontogeny
# (island_ontogeny = 0, sea_level = 1).
#
# The simulation (island_area()) uses an area that oscillates with the sea
# level. The likelihood used to ignore it: island_area_vector() returned area 1
# as soon as island_ontogeny was 0, before the sea level was looked at (same
# shortcut in the C++ island_area_vector() in src/DAISIE_time_rates.h), and
# DAISIE_loglik_CS_M1() scaled K with the maximum area pars1[8] (which is not the
# maximum area of such an island: its largest area is current_area = pars1[9]).
# DAISIE_ML accepts island_ontogeny = 0 with sea_level = 1, so ML on
# sea-level-only models ignored the sea level.

sl_area_pars <- create_area_pars(
  max_area = 1, current_area = 500, proportional_peak_t = 0.5,
  total_island_age = 5, sea_level_amplitude = 3, sea_level_frequency = 2,
  island_gradient_angle = 45
)
age <- 5
ap_vec <- function(ap) {
  unlist(ap)[c("max_area", "current_area", "proportional_peak_t",
               "total_island_age", "sea_level_amplitude",
               "sea_level_frequency", "island_gradient_angle")]
}
# loglik of one clade of a sea-level-only island
ll_clade_sea <- function(p5, d, x, ap, brts, stac, miss, methode) {
  DAISIE_loglik_CS_M1(
    pars1 = c(p5, d, x, ap_vec(ap), 0, 1, age, 1),
    pars2 = c(30, 11, 0, 0, 0, 1), brts = brts, stac = stac,
    missnumspec = miss, methode = methode
  )
}
clade_cases <- list(
  list(stac = 0, brts = c(age, 0),         miss = 0),
  list(stac = 2, brts = c(age, 3, 0),      miss = 0),
  list(stac = 3, brts = c(age, 4, 2, 0),   miss = 2),
  list(stac = 4, brts = c(age, 3, 0),      miss = 0),
  list(stac = 1, brts = c(age, age, 0),    miss = 3),
  list(stac = 5, brts = c(age, 3.5, 0),    miss = 2),
  list(stac = 6, brts = c(age, 4, 2, 0),   miss = 2),
  list(stac = 9, brts = c(age, 4, 2, 0),   miss = 2)
)

test_that("island_area_vector() agrees with island_area() for sea level only", {
  for (tt in c(0, 1, 2.5, 4, 5)) {
    expect_equal(
      suppressWarnings(island_area_vector(
        timeval = tt, area_pars = unlist(sl_area_pars)[1:7],
        island_ontogeny = 0, sea_level = 1, total_time = 5, peak = 1
      )),
      island_area(
        timeval = tt, total_time = 5, area_pars = sl_area_pars, peak = 1,
        island_ontogeny = 0, sea_level = 1
      ),
      info = paste("t =", tt)
    )
  }
})

test_that("sea level with amplitude 0 gives the time-constant model with constant area", {
  # amplitude 0: the area is constant, equal to current_area A0, so the
  # time-dependent model is the time-constant model with rates
  # lac * A0^d, mu * A0^(-x) and carrying capacity K * A0.
  for (A0 in c(7, 40)) {
    ap0 <- create_area_pars(
      max_area = 1, current_area = A0, proportional_peak_t = 0.5,
      total_island_age = age, sea_level_amplitude = 0, sea_level_frequency = 2,
      island_gradient_angle = 45
    )
    p5 <- c(0.5, 0.3, 4, 0.1, 0.4)
    d <- 0.3
    x <- 0.4
    p5_const <- c(p5[1] * A0^d, p5[2] * A0^(-x), p5[3] * A0, p5[4], p5[5])
    for (cs in clade_cases) {
      expect_equal(
        ll_clade_sea(p5, d, x, ap0, cs$brts, cs$stac, cs$miss,
                     "odeint::runge_kutta_dopri5"),
        DAISIE_loglik_CS_M1(
          pars1 = p5_const, pars2 = c(30, 11, 0, 0, NA, NA), brts = cs$brts,
          stac = cs$stac, missnumspec = cs$miss,
          methode = "odeint::runge_kutta_dopri5"
        ),
        tolerance = 1e-10,
        info = paste("A0 =", A0, "stac =", cs$stac, "miss =", cs$miss)
      )
    }
  }
})

test_that("R (lsodes) and C++ (odeint) loglik agree for an oscillating sea level", {
  skip_on_cran()
  p5 <- c(0.5, 0.3, 0.05, 0.1, 0.4)
  for (cs in clade_cases[1:4]) {
    expect_equal(
      ll_clade_sea(p5, 0.3, 0.3, sl_area_pars, cs$brts, cs$stac, cs$miss,
                   "lsodes"),
      ll_clade_sea(p5, 0.3, 0.3, sl_area_pars, cs$brts, cs$stac, cs$miss,
                   "odeint::runge_kutta_dopri5"),
      tolerance = 1e-8,
      info = paste("stac =", cs$stac, "miss =", cs$miss)
    )
  }
})

test_that("loglik depends on sea level for a sea-level-only island", {
  skip_on_cran()
  set.seed(2)
  dl <- DAISIE_sim_cr(time = 5, M = 40, pars = c(0.5, 0.3, 20, 0.1, 0.4),
                      replicates = 1, plot_sims = FALSE, verbose = FALSE)[[1]]
  ap <- unlist(sl_area_pars)[1:7]
  ll <- function(sea_level, amplitude) {
    ap[5] <- amplitude
    suppressWarnings(DAISIE_loglik_all(
      pars1 = c(0.5, 0.3, 20, 0.1, 0.4, 0.3, 0.3, ap),
      pars2 = c(40, 11, 0, 0, 0, sea_level), datalist = dl,
      methode = "odeint::runge_kutta_dopri5"
    ))
  }
  expect_true(abs(ll(sea_level = 1, amplitude = 3) -
                    ll(sea_level = 0, amplitude = 3)) > 1e-3)
  expect_true(abs(ll(sea_level = 1, amplitude = 3) -
                    ll(sea_level = 1, amplitude = 0)) > 1e-3)
})

test_that("loglik of an empty island matches the simulation for a sea-level-only island", {
  skip_on_cran()
  # d = 0, x = 0.5: extinction is strongly area dependent, so the probability
  # that the island is empty differs clearly between an island with sea level
  # change (area oscillating between 137 and 500) and one with a constant area
  # of 500 (0.35 vs 0.28 here).
  p5 <- c(0.5, 6, 0.05, 0.4, 0.5)
  reps <- 600
  set.seed(11)
  sim <- DAISIE_sim_time_dep(
    time = age, M = 1, pars = p5, replicates = reps, island_ontogeny = "const",
    sea_level = "sine", area_pars = sl_area_pars,
    hyper_pars = create_hyper_pars(d = 0, x = 0.5), plot_sims = FALSE,
    verbose = FALSE
  )
  frac_empty <- mean(vapply(sim, function(r) r[[1]]$not_present == 1, logical(1)))
  p_empty <- exp(ll_clade_sea(p5, 0, 0.5, sl_area_pars, c(age, 0), 0, 0,
                              "odeint::runge_kutta_dopri5"))
  expect_lt(abs(frac_empty - p_empty), 4 * sqrt(p_empty * (1 - p_empty) / reps))
})
