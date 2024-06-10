# Tests for the time-dependent likelihood computation

# Mock arguments shared across models
initprobs <- rep(0.2, 5)
tvec <- c(-2, -1)
rtol <- 1E-10
atol <- 1E-16
method <- "deSolve_R::ode45"

# Pick some parameter values for the basic model
pars <- create_pars(
  clado_rate = 0.1,
  ext_rate = 0,
  carr_cap = 100,
  immig_rate = 0.01,
  ana_rate = 0.01
)

# Extra parameters relevant for the function being tested
extra_pars <- list(kk = 3, ddep = 11)

# Island age
total_time <- 4

# Note: ddep = 11 is to make sure the time-independent version is run with
# diversity-dependent immigration rates (which is the default in the time-
# dependent algorithm, see function get_immig_rate_per_capita()).

# Mock area parameters
area_pars <- create_area_pars(
  max_area = 1,
  current_area = 0.9,
  proportional_peak_t = 0,
  total_island_age = 5, # island lifespan
  sea_level_amplitude = 0,
  sea_level_frequency = 0,
  island_gradient_angle = 0
)

# Assemble parameters for the time-independent model
pars_ti <- c(pars, extra_pars)

# And for the time-dependent model
pars_td <- c(
  pars,
  create_hyper_pars(d = 1, x = 1),
  ka = 1,
  area_pars,
  island_ontogeny = 1, # ontogeny is on
  sea_level = 0,
  total_time = total_time,
  peak = DAISIE:::calc_peak(total_time, area_pars),
  extra_pars
)

# Unpack both into vectors and not lists
pars_ti <- unlist(pars_ti)
pars_td <- unlist(pars_td)

# Time-dependent and independent are the same with no island ontogeny
test_that("time-dependent integration with no ontogeny", {

  # Make sure ontogeny is off
  pars_td["island_ontogeny"] <- 0

  # Compute time-independent likelihood
  loglik_ti <- DAISIE:::DAISIE_integrate_const(
    initprobs = initprobs, tvec = tvec, rhs_func = DAISIE:::DAISIE_loglik_rhs,
    pars = pars_ti, rtol = rtol, atol = atol, method = "deSolve_R::ode45"
  )

  # Compute time-dependent likelihood
  loglik_td <- DAISIE:::DAISIE_integrate_time(
    initprobs = initprobs, tvec = tvec, rhs_func = DAISIE:::DAISIE_loglik_rhs,
    pars = pars_td, rtol = rtol, atol = atol, method = "ode45"
  )

  # Check that both are equal
  expect_true(all(loglik_td == loglik_ti))

})

# Both are the same with ontogeny but constant area
test_that("time-dependent integration with constant area", {

  # Note: ontogeny should be on

  # Make sure the island area is constant through time
  pars_td["peak"] <- 0

  # Compute time-independent likelihood
  loglik_ti <- DAISIE:::DAISIE_integrate_const(
    initprobs = initprobs, tvec = tvec, rhs_func = DAISIE:::DAISIE_loglik_rhs,
    pars = pars_ti, rtol = rtol, atol = atol, method = "deSolve_R::ode45"
  )

  # Compute time-dependent likelihood
  loglik_td <- DAISIE:::DAISIE_integrate_time(
    initprobs = initprobs, tvec = tvec, rhs_func = DAISIE:::DAISIE_loglik_rhs,
    pars = pars_td, rtol = rtol, atol = atol, method = "ode45"
  )

  # Check that both are equal
  expect_true(all(loglik_td == loglik_ti))

})

# Both are almost the same when island area changes very little through time
test_that("time-dependent integration with slight island ontogeny", {

  # Note: ontogeny should be on

  # Compute time-independent likelihood
  loglik_ti <- DAISIE:::DAISIE_integrate_const(
    initprobs = initprobs, tvec = tvec, rhs_func = DAISIE:::DAISIE_loglik_rhs,
    pars = pars_ti, rtol = rtol, atol = atol, method = "deSolve_R::ode45"
  )

  # Compute time-dependent likelihood
  loglik_td <- DAISIE:::DAISIE_integrate_time(
    initprobs = initprobs, tvec = tvec, rhs_func = DAISIE:::DAISIE_loglik_rhs,
    pars = pars_td, rtol = rtol, atol = atol, method = "ode45"
  )

  # Compute some error metric
  errors <- (loglik_td - loglik_ti)^2

  # Check that it is small
  expect_true(all(errors < 1e-5))

})

# Both should diverge more when ontogeny becomes more prominent
test_that("time-dependent integration with increasingly important ontogeny", {

  # Note: we are testing for a "dose-response relationship" as indication
  # that time-dependence has been correctly implemented.

  # Prepare three islands with increasingly weaker peakiness
  pars_td1 <- pars_td2 <- pars_td
  pars_td1["peak"] <- pars_td["peak"] / 2
  pars_td2["peak"] <- pars_td1["peak"] / 2

  # Note: the smaller the peakiness the weaker the ontogenic change through
  # time.

  # Compute time-independent likelihood
  loglik_ti <- DAISIE:::DAISIE_integrate_const(
    initprobs = initprobs, tvec = tvec, rhs_func = DAISIE:::DAISIE_loglik_rhs,
    pars = pars_ti, rtol = rtol, atol = atol, method = "deSolve_R::ode45"
  )

  # Compute time-dependent likelihood
  loglik_td <- DAISIE:::DAISIE_integrate_time(
    initprobs = initprobs, tvec = tvec, rhs_func = DAISIE:::DAISIE_loglik_rhs,
    pars = pars_td, rtol = rtol, atol = atol, method = "ode45"
  )

  # And again...
  loglik_td1 <- DAISIE:::DAISIE_integrate_time(
    initprobs = initprobs, tvec = tvec, rhs_func = DAISIE:::DAISIE_loglik_rhs,
    pars = pars_td1, rtol = rtol, atol = atol, method = "ode45"
  )

  # ... and again
  loglik_td2 <- DAISIE:::DAISIE_integrate_time(
    initprobs = initprobs, tvec = tvec, rhs_func = DAISIE:::DAISIE_loglik_rhs,
    pars = pars_td2, rtol = rtol, atol = atol, method = "ode45"
  )

  # Compute errors
  errors <- (loglik_td - loglik_ti)^2
  errors1 <- (loglik_td1 - loglik_ti)^2
  errors2 <- (loglik_td2 - loglik_ti)^2

  # Check that deviations are smaller when peakiness is weaker
  expect_true(all(errors >= errors1))
  expect_true(all(errors1 >= errors2))

})
