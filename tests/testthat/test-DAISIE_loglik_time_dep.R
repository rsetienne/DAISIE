# ---------------------------------------------------------------------------
# DAISIE_loglik (constant rates) compared with DAISIE_loglik_time
# (time-dependent rates) with very low island ontogeny and sea level effects.
# ---------------------------------------------------------------------------

library(DAISIE)

test_that("island ontogeny and sea level with low values of d and x: DAISIE_loglik_CS
          equals the constant-rate loglik", {
            lac <- 0.5
            mu  <- 0.3
            K   <- Inf
            gam <- 0.02
            laa <- 1.0
            res <- 100

            d <- 0.000001
            x <- 0.000001
            island_ontogeny <- 1
            sea_level <- 1

            utils::data("Galapagos_datalist", package = "DAISIE", envir = environment())
            datalist <- Galapagos_datalist
            island_age <- datalist[[1]]$island_age

            area_pars <- create_area_pars(
              max_area              = 1000,
              current_area          = 500,
              proportional_peak_t   = 0.2,
              total_island_age      = 2 * island_age,
              sea_level_amplitude   = 1,
              sea_level_frequency   = 2,
              island_gradient_angle = 45
            )

            pars1_time <- c(lac, mu, K, gam, laa,
                            d, x,
                            area_pars$max_area,
                            area_pars$current_area,
                            area_pars$proportional_peak_t,
                            area_pars$total_island_age,
                            area_pars$sea_level_amplitude,
                            area_pars$sea_level_frequency,
                            area_pars$island_gradient_angle)

            loglik_cr <- DAISIE_loglik_CS(
              pars1 = c(lac, mu, K, gam, laa),
              pars2 = c(res, 11, 0, 0),
              datalist = datalist,
              methode = "lsodes",
              CS_version = list(model = 1, function_to_optimize = 'DAISIE', sampling = 'n'),
              abstolint = 1E-16,
              reltolint = 1E-10)

            loglik_time <- DAISIE_loglik_CS(
              pars1 = pars1_time,
              pars2 = c(res, 11, 0, 0, island_ontogeny, sea_level),
              datalist = datalist,
              methode = "lsodes",
              CS_version = list(model = 1, function_to_optimize = 'DAISIE', sampling = 'n'),
              abstolint = 1E-16,
              reltolint = 1E-10)

            testthat::expect_lt(abs(loglik_time - loglik_cr), 1e-3)
          })
