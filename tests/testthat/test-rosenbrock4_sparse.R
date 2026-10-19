# The 'odeint::rosenbrock4' stepper uses the sparse (banded) Rosenbrock method
# of src/DAISIE_stiff.h. deSolve's lsodes serves as reference.

cs_init <- function(n, seed = 3) {
  set.seed(seed)
  x <- runif(n)
  x / sum(x)
}

test_that("sparse rosenbrock4 equals lsodes for all CS systems", {
  lx <- 15
  pars_model <- c(0.5, 0.3, 8, 0.05, 0.4)
  for (kk in 0:1) {
    for (rm in list(
      list(rhs = DAISIE_loglik_rhs, n = 2 * lx + 1),
      list(rhs = DAISIE_loglik_rhs1, n = 4 * lx),
      list(rhs = DAISIE_loglik_rhs2, n = 3 * lx)
    )) {
      init <- cs_init(rm$n)
      pars <- c(pars_model, kk, 11)
      y_ref <- DAISIE:::DAISIE_integrate(
        init, c(-4, 0), rm$rhs, pars, 1e-13, 1e-13, "lsodes"
      )
      y_new <- DAISIE:::DAISIE_integrate(
        init, c(-4, 0), rm$rhs, pars, 1e-12, 1e-14, "odeint::rosenbrock4"
      )
      testthat::expect_equal(unname(y_new), unname(y_ref), tolerance = 1e-6,
                             info = paste("kk =", kk, "n =", rm$n))
    }
  }
})

test_that("sparse rosenbrock4 equals deSolve on a stiff CS problem", {
  lx <- 40
  pars <- c(30, 2.68, 10, 0.0093, 1.01, 0, 11)   # high cladogenesis: stiff
  init <- c(1, rep(0, 2 * lx))
  y_r <- DAISIE:::DAISIE_integrate(
    init, c(-4, 0), DAISIE_loglik_rhs, pars, 1e-14, 1e-12, "lsodes"
  )
  y_new <- DAISIE:::DAISIE_integrate(
    init, c(-4, 0), DAISIE_loglik_rhs, pars, 1e-14, 1e-12, "odeint::rosenbrock4"
  )
  testthat::expect_equal(unname(y_new), unname(y_r), tolerance = 1e-6)
})

test_that("rosenbrock4 loglik equals lsodes, two types", {
  utils::data(Galapagos_datalist_2types)
  pars1 <- c(0.195442017, 0.087959583, Inf, 0.002247364, 0.873605049,
             3755.202241, 8.909285094, 14.99999923, 0.002247364,
             0.873605049, 0.163)
  pars2 <- c(40, 11, 0, 1)
  loglik_lsodes <- DAISIE_loglik_all(
    pars1, pars2, Galapagos_datalist_2types, methode = "lsodes"
  )
  loglik_rb <- DAISIE_loglik_all(
    pars1, pars2, Galapagos_datalist_2types, methode = "odeint::rosenbrock4"
  )
  testthat::expect_equal(loglik_rb, loglik_lsodes, tolerance = 1e-7)
})

test_that("rosenbrock4 works for DAISIE_DE", {
  utils::data(made_up_datalist)
  pars1 <- c(2.546591, 2.678781, 2.678781, 0.009326754, 1.008583)
  ll <- function(methode) {
    DAISIE:::DAISIE_loglik_CS(
      pars1 = pars1, pars2 = c(100, 11, 0, 0), datalist = made_up_datalist,
      methode = methode, abstolint = 1e-14, reltolint = 1e-12,
      CS_version = list(model = 1, function_to_optimize = "DAISIE_DE",
                        sampling = "n")
    )
  }
  testthat::expect_equal(ll("odeint::rosenbrock4"),
                         ll("odeint::runge_kutta_cash_karp54"),
                         tolerance = 1e-8)
})

test_that("rosenbrock4 works for the IW model (sparse Jacobian and LU)", {
  skip_on_cran()
  utils::data(frogs_datalist, package = "DAISIE")
  pars1 <- c(0.2, 0.1, 1000.1, 0.001, 0.3)
  pars2 <- c(10, 11, 0, 0)
  ll <- function(methode) {
    DAISIE_loglik_IW(pars1 = pars1, pars2 = pars2, datalist = frogs_datalist,
                     methode = methode, abstolint = 1e-12, reltolint = 1e-10)
  }
  testthat::expect_equal(ll("odeint::rosenbrock4"),
                         ll("odeint::runge_kutta_cash_karp54"),
                         tolerance = 1e-7)
})
