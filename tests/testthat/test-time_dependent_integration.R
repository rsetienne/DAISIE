# Tests for the time-dependent likelihood computation

# Time-dependent and independent are the same with no island ontogeny
test_that("time-dependent integration with no ontogeny", {

  # Hands-on tutorial would work I think. Perhaps until then get familiar with
  # the integration procedure. That is, both the paper and the extra equations.

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

  # Note: ddep = 11 is to make sure the time-independent version is run with
  # diversity-dependent immigration rates (which is the default in the time-
  # dependent algorithm, see function get_immig_rate_per_capita()).

  # Island age
  total_time <- 4

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
    ka = 2,
    area_pars,
    island_ontogeny = 0,
    sea_level = 0,
    total_time = total_time,
    peak = DAISIE:::calc_peak(total_time, area_pars),
    extra_pars
  )

  # Unpack both into vectors and not lists
  pars_ti <- unlist(pars_ti)
  pars_td <- unlist(pars_td)

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
  expect_equal(loglik_td, loglik_ti)

})
