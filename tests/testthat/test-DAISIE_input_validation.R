test_that("calc_next_timeval_shift returns a numeric dt when all rates are zero", {
  out <- DAISIE:::calc_next_timeval_shift(
    max_rates = list(0, 0, 0, 0),
    timeval = 1,
    dynamic_shift_times = c(3, 4),
    total_time = 5
  )
  testthat::expect_true(is.numeric(out$dt))
  testthat::expect_equal(out$dt, 4)
  testthat::expect_equal(out$timeval, 3) # shift time comes first
  testthat::expect_true(out$rate_shift)
})

test_that("branching times older than the island age are rejected", {
  age <- 5
  pars <- c(0.5, 0.3, 15, 0.1, 0.4)
  ll <- function(brts) {
    DAISIE_loglik_all(
      pars1 = pars, pars2 = c(30, 11, 0, 0),
      datalist = list(
        list(island_age = age, not_present = 1),
        list(colonist_name = "a", branching_times = brts, stac = 2,
             missing_species = 0, type1or2 = 1)
      ),
      methode = "lsodes"
    )
  }
  testthat::expect_true(is.finite(ll(c(age, 3))))
  # colonisation exactly at the island age is allowed
  testthat::expect_true(is.finite(ll(c(age, age))))
  testthat::expect_error(ll(c(age, age + 1)), "older than the island age")
})

test_that("beta ontogeny needs an island age below the total island age", {
  age <- 5
  pars <- c(0.5, 0.3, 15, 0.1, 0.4, 0.1, 0.1, 1000, 500, 0.1, age, 0, 0, 0)
  dl <- list(
    list(island_age = age, not_present = 1),
    list(colonist_name = "a", branching_times = c(age, 3), stac = 2,
         missing_species = 0, type1or2 = 1)
  )
  testthat::expect_error(
    DAISIE_loglik_all(pars1 = pars, pars2 = c(30, 11, 0, 0, 1, 0),
                      datalist = dl, methode = "lsodes"),
    "must be smaller than the total island age"
  )
})

test_that("constant ontogeny does not need a peak, also at total island age", {
  age <- 5
  pars_static <- c(0.5, 0.3, 15, 0.1, 0.4)
  dl <- list(
    list(island_age = age, not_present = 1),
    list(colonist_name = "a", branching_times = c(age, 3), stac = 2,
         missing_species = 0, type1or2 = 1)
  )
  ll_static <- DAISIE_loglik_all(
    pars1 = pars_static, pars2 = c(30, 11, 0, 0), datalist = dl,
    methode = "lsodes"
  )
  # total island age == island age, constant area
  pars_time <- c(pars_static, 0, 0, 1, 1, 0, age, 0, 0, 0)
  ll_time <- DAISIE_loglik_all(
    pars1 = pars_time, pars2 = c(30, 11, 0, 0, 0, 0), datalist = dl,
    methode = "lsodes"
  )
  testthat::expect_equal(ll_time, ll_static, tolerance = 1e-6)
})
