# Tests of the C++ (boost::odeint) integrator of the time-dependent CS model,
# daisie_odeint_cs_time(), against the R/deSolve implementation.

cpp_time_pars <- function(sea_level = 0, amplitude = 0, frequency = 0,
                          island_age = 5, total_age = 10) {
  area_pars <- create_area_pars(1000, 500, 0.1, total_age, amplitude,
                                frequency, if (sea_level == 1) 10 else 0)
  peak <- calc_peak(island_age, area_pars)
  c(0.5, 0.3, 20, 0.15, 0.8, 0.1, 0.1,
    1000, 500, 0.1, total_age, amplitude, frequency, if (sea_level == 1) 10 else 0,
    1, sea_level, island_age, peak)
}

cpp_time_runmods <- function(lx) {
  list(
    list(rhs = DAISIE_loglik_rhs, n = 2 * lx + 1),
    list(rhs = DAISIE_loglik_rhs1, n = 4 * lx),
    list(rhs = DAISIE_loglik_rhs2, n = 3 * lx)
  )
}

test_that("C++ time-dependent integration equals deSolve for all rhs and kk", {
  lx <- 10
  pars_model <- cpp_time_pars()
  set.seed(7)
  for (kk in 0:1) {
    for (rm in cpp_time_runmods(lx)) {
      init <- runif(rm$n)
      init <- init / sum(init)
      pars <- c(pars_model, kk, 11)
      y_r <- DAISIE:::DAISIE_integrate(
        init, c(-5, 0), rm$rhs, pars, 1e-12, 1e-16, "lsodes"
      )
      print(y_r)
      for (stepper in c("odeint::runge_kutta_dopri5",
                        "odeint::runge_kutta_fehlberg78",
                        "odeint::runge_kutta_cash_karp54",
                        "odeint::bulirsch_stoer")) {
        y_cpp <- DAISIE:::DAISIE_integrate(
          init, c(-5, 0), rm$rhs, pars, 1e-12, 1e-16, stepper
        )
        testthat::expect_equal(
          unname(y_cpp), unname(y_r), tolerance = 1e-8,
          info = paste("kk =", kk, "n =", rm$n, stepper)
        )
      }
    }
  }
})

test_that("C++ time-dependent integration handles ontogeny with sea level", {
  lx <- 10
  pars <- c(cpp_time_pars(sea_level = 1, amplitude = 1, frequency = 2), 0, 11)
  init <- c(1, rep(0, 2 * lx))
  y_r <- DAISIE:::DAISIE_integrate(
    init, c(-5, 0), DAISIE_loglik_rhs, pars, 1e-12, 1e-16, "lsodes"
  )
  y_cpp <- DAISIE:::DAISIE_integrate(
    init, c(-5, 0), DAISIE_loglik_rhs, pars, 1e-12, 1e-16,
    "odeint::runge_kutta_dopri5"
  )
  testthat::expect_equal(unname(y_cpp), unname(y_r), tolerance = 1e-8)
  # the sea level changes the solution
  pars_no_sea <- c(cpp_time_pars(), 0, 11)
  y_no_sea <- DAISIE:::DAISIE_integrate(
    init, c(-5, 0), DAISIE_loglik_rhs, pars_no_sea, 1e-12, 1e-16,
    "odeint::runge_kutta_dopri5"
  )
  testthat::expect_gt(max(abs(y_cpp - y_no_sea)), 1e-6)
})

test_that("C++ time-dependent integration reduces to the C++ constant-area one", {
  lx <- 10
  pars_static <- c(0.5, 0.3, 20, 0.05, 0.8)
  pars_time <- c(pars_static, 0, 0, 1, 1, 0, 10, 0, 0, 0, 0, 0, 10, 1)
  set.seed(3)
  for (rm in cpp_time_runmods(lx)) {
    init <- runif(rm$n)
    init <- init / sum(init)
    y_const <- DAISIE:::DAISIE_integrate(
      init, c(-10, 0), rm$rhs, c(pars_static, 0, 11), 1e-12, 1e-16,
      "odeint::runge_kutta_dopri5"
    )
    y_time <- DAISIE:::DAISIE_integrate(
      init, c(-10, 0), rm$rhs, c(pars_time, 0, 11), 1e-12, 1e-16,
      "odeint::runge_kutta_dopri5"
    )
    testthat::expect_equal(unname(y_time), unname(y_const), tolerance = 1e-9,
                           info = paste("n =", rm$n))
  }
})

test_that("rosenbrock4 works for the time-dependent model", {
  lx <- 6
  pars <- c(cpp_time_pars(), 0, 11)
  init <- c(1, rep(0, 2 * lx))
  y_r <- DAISIE:::DAISIE_integrate(
    init, c(-5, 0), DAISIE_loglik_rhs, pars, 1e-12, 1e-16, "lsodes"
  )
  y_ros <- DAISIE:::DAISIE_integrate(
    init, c(-5, 0), DAISIE_loglik_rhs, pars, 1e-6, 1e-9,
    "odeint::rosenbrock4"
  )
  testthat::expect_equal(unname(y_ros), unname(y_r), tolerance = 1e-4)
})

test_that("time-dependent likelihood is the same with deSolve and odeint", {
  pars <- cpp_time_pars()
  age <- 5
  cases <- list(
    list(stac = 0, brts = c(age, 0)),
    list(stac = 2, brts = c(age, 3, 0)),
    list(stac = 3, brts = c(age, 4, 2, 0)),
    list(stac = 1, brts = c(age, age, 0)),
    list(stac = 5, brts = c(age, age, 0)),
    list(stac = 6, brts = c(age, age, 2, 0)),
    list(stac = 9, brts = c(age, 3, 2, 0))
  )
  for (case in cases) {
    ll <- function(methode) {
      DAISIE_loglik_CS_M1(
        pars1 = pars, pars2 = c(30, 11, 0, 0, 1, 0), brts = case$brts,
        stac = case$stac, missnumspec = 1, methode = methode,
        reltolint = 1e-12, abstolint = 1e-16
      )
    }
    testthat::expect_equal(
      ll("odeint::runge_kutta_dopri5"), ll("lsodes"), tolerance = 1e-7,
      info = paste("stac =", case$stac)
    )
  }
})

test_that("daisie_odeint_cs_time rejects invalid input", {
  init <- c(1, rep(0, 8))
  testthat::expect_error(
    DAISIE_ode_cs_time(init, c(-5, 0), c(cpp_time_pars(), 0, 11), 1e-10, 1e-10,
                       "odeint::runge_kutta_dopri5", runmod = "nonsense"),
    "unknown runmod"
  )
  testthat::expect_error(
    .Call("daisie_odeint_cs_time", "daisie_runmod", init, c(-5, 0), 4, 0,
          cpp_time_pars()[1:10], "odeint::runge_kutta_dopri5", 1e-10, 1e-10),
    "18 model parameters"
  )
})
