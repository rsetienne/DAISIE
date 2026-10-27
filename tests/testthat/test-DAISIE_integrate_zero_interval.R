# A colonisation at the island age (e.g. a maximum colonisation time equal to
# the island age) gives an integration interval of length 0, e.g. c(-10, -10).
# deSolve::ode(method = "ode45") must not be called with such an interval: it
# then uses hmax = .Machine$double.xmax, evaluates the right-hand side at times
# of the order of 1e307 and lets the state overflow to Inf/NaN unless the
# derivative happens to be exactly 0 (as it was for beta-function ontogeny
# only because the area is 0 at the moment of island emergence and the
# resulting NaN was turned into a rate of 0). The solution on an interval of
# length 0 is the initial state.

island_age <- 10
pars_beta <- c(
  0.5, 0.3, 20, 0.05, 0.8, 0.1, 0.1, # lac, mu, K, gam, laa, d, x
  1, 0.5, 0.5, island_age, 0, 0, 0,  # area pars
  1, 0, island_age, 0.5              # ontogeny, sea level, total time, peak
)

test_that("DAISIE_integrate returns the initial state on an interval of length 0", {
  lx <- 12
  init <- c(0.3, rep(0.05, 2 * lx))
  for (method in c("ode45", "lsodes", "lsoda")) {
    expect_equal(
      DAISIE_integrate(
        init, c(-island_age, -island_age), DAISIE_loglik_rhs,
        c(pars_beta, 0, 11), rtol = 1e-10, atol = 1e-14, method = method
      ),
      init,
      info = paste("time-dependent,", method)
    )
    expect_equal(
      DAISIE_integrate(
        init, c(-island_age, -island_age), DAISIE_loglik_rhs,
        c(0.5, 0.3, 20, 0.05, 0.8, 0, 11), rtol = 1e-10, atol = 1e-14,
        method = method
      ),
      init,
      info = paste("time-constant,", method)
    )
  }
})

test_that("ode45 gives a finite loglik when the colonisation is at the island age", {
  skip_on_cran()
  ll <- function(stac, brts, methode) {
    DAISIE_loglik_CS_M1(
      pars1 = pars_beta, pars2 = c(30, 11, 0, 0, 1, 0), brts = brts,
      stac = stac, missnumspec = 0, methode = methode,
      reltolint = 1e-10, abstolint = 1e-14
    )
  }
  for (case in list(
    list(stac = 1, brts = c(island_age, island_age, 0)),
    list(stac = 5, brts = c(island_age, island_age, 0)),
    list(stac = 6, brts = c(island_age, island_age, 3, 0))
  )) {
    ref <- ll(case$stac, case$brts, "lsodes")
    expect_true(is.finite(ref))
    expect_equal(ll(case$stac, case$brts, "ode45"), ref, tolerance = 1e-6,
                 info = paste("stac", case$stac))
  }
})
