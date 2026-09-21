# ---------------------------------------------------------------------------
# Time-dependent DAISIE compared with constant-rate DAISIE on the Galapagos
# data, with island ontogeny and sea level.
# DAISIE_loglik_CS is time-dependent when it gets 14 parameters and
# constant-rate when it gets 5.
#
# Use methode = "lsodes": the time-dependent likelihood does not support odeint yet.
# ---------------------------------------------------------------------------

lac <- 0.5      # cladogenesis
mu  <- 0.3      # extinction
K   <- Inf      # carrying capacity
gam <- 0.02     # immigration
laa <- 1.0      # anagenesis
res <- 30

utils::data(Galapagos_datalist, package = "DAISIE", envir = environment())
island_age <- Galapagos_datalist[[1]]$island_age   # 4

# constant-rate loglik of the Galapagos data
loglik_cr <- DAISIE_loglik_CS(
  pars1 = c(lac, mu, K, gam, laa),
  pars2 = c(res, 11, 0, 0),                # res, ddmodel, cond, verbose
  datalist = Galapagos_datalist,
  methode = "lsodes"
)

test_that("constant area: time-dependent loglik equals constant-rate loglik", {
  # the area is always 1, so d and x have no effect
  loglik_time <- DAISIE_loglik_CS(
    pars1 = c(lac, mu, K, gam, laa,
              0.007, 0.009,                # d, x
              1000, 500, 0.2, 10,          # max_area, current_area, proportional_peak_t, total_island_age
              1, 3, 45),                    # sea_level_amplitude, sea_level_frequency, island_gradient_angle
    pars2 = c(res, 11, 0, 0, 0, 0),        # res, ddmodel, cond, verbose, island_ontogeny, sea_level
    datalist = Galapagos_datalist,
    methode = "lsodes"
  )
  # total_island_age = 2 * island_age puts the present at the area peak;
  # otherwise calc_peak() fails when max_area = current_area
  testthat::expect_lt(abs(loglik_time - loglik_cr), 1e-7)
})

test_that("island ontogeny with d = 0 and x = 0: time-dependent loglik equals
          constant-rate loglik", {
            loglik_time <- DAISIE_loglik_CS(
              pars1 = c(lac, mu, K, gam, laa,
                        0, 0,                        # d, x
                        1000, 500, 0.2, 10,          # max_area, current_area, proportional_peak_t, total_island_age
                        0, 0, 0),                    # no sea level
              pars2 = c(res, 11, 0, 0, 1, 0),        # island_ontogeny = 1, sea_level = 0
              datalist = Galapagos_datalist,
              methode = "lsodes"
            )
            testthat::expect_lt(abs(loglik_time - loglik_cr), 1e-7)
          })

test_that("island ontogeny and sea level with d = 0 and x = 0: logliks are close", {
  # Not exactly equal: right after the island emerges its area is 0, and at
  # area 0 there is no immigration. Each of the 992 mainland species that did
  # not colonise adds a little, so the island loglik differs by about 0.1.
  loglik_time <- DAISIE_loglik_CS(
    pars1 = c(lac, mu, K, gam, laa,
              0, 0,                        # d, x
              1000, 500, 0.2, 10,          # max_area, current_area, proportional_peak_t, total_island_age
              1, 3, 45),                   # sea_level_amplitude, sea_level_frequency, island_gradient_angle
    pars2 = c(res, 11, 0, 0, 1, 1),        # island_ontogeny = 1, sea_level = 1
    datalist = Galapagos_datalist,
    methode = "lsodes"
  )
  testthat::expect_lt(abs(loglik_time - loglik_cr), 1)
})

test_that("island ontogeny and sea level with d and x: loglik is finite and
          differs from constant-rate loglik", {
            loglik_time <- DAISIE_loglik_CS(
              pars1 = c(lac, mu, K, gam, laa,
                        0.007, 0.009,                # d, x
                        1000, 500, 0.2, 10,          # max_area, current_area, proportional_peak_t, total_island_age
                        1, 3, 45),                   # sea_level_amplitude, sea_level_frequency, island_gradient_angle
              pars2 = c(res, 11, 0, 0, 1, 1),        # island_ontogeny = 1, sea_level = 1
              datalist = Galapagos_datalist,
              methode = "lsodes"
            )
            testthat::expect_true(is.finite(loglik_time))
            testthat::expect_gt(abs(loglik_time - loglik_cr), 1e-3)
          })
