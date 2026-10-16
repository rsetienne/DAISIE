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
