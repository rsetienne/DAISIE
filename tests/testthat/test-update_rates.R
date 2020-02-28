context("update_rates")

test_that("update_rates constant rates is silent and gives correct output", {
  timeval <- 0
  totaltime <- 1
  gam <- 0.009
  mu <- 2.0
  laa <- 1.0
  lac <- 2.5
  K <- 3
  default_pars <- create_default_pars(
    island_ontogeny = 0,
    sea_level = 0,
    area_pars = NULL,
    hyper_pars = NULL,
    dist_pars = NULL,
    ext_pars = NULL,
    totaltime = totaltime,
    pars = c(lac, mu, K, gam, laa)
  )
  island_ontogeny <- translate_island_ontogeny("const")
  sea_level <- translate_sea_level("const")
  extcutoff <- 1000.0
  num_spec <- 0
  num_immigrants <- 0
  mainland_n <- 1
  set.seed(42)
  expect_silent(rates <- update_rates(
    timeval = timeval,
    totaltime = totaltime,
    gam = gam,
    laa = laa,
    lac = lac,
    mu = mu,
    hyper_pars = default_pars$hyper_pars,
    area_pars = default_pars$area_pars,
    dist_pars = default_pars$dist_pars,
    ext_pars = default_pars$ext_pars,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    extcutoff = extcutoff,
    K = K,
    num_spec = num_spec,
    num_immigrants = num_immigrants,
    mainland_n = mainland_n))
  expect_true(are_rates(rates))
  expected_rates <- list(immig_rate = 0.009,
                         ext_rate = 0,
                         ana_rate = 0,
                         clado_rate = 0)
  expect_equal(rates, expected_rates)
})


test_that("update area-dependent rates is silent and gives correct output", {
  set.seed(42)
  lac <- 2.5
  laa <- 1.0
  gam <- 0.009
  mu <- 2.5
  K <- 3
  default_pars <- create_default_pars(
    island_ontogeny = 1,
    sea_level = 0,
    area_pars = create_area_pars(
      max_area = 1.0,
      proportional_peak_t = 0.5,
      peak_sharpness = 1.0,
      total_island_age = 1.0,
      sea_level_amplitude = 0,
      sea_level_frequency = 0,
      island_gradient_angle = 0),
    hyper_pars = NULL,
    dist_pars = NULL,
    ext_pars = c(0.5, 10.0),
    totaltime = totaltime,
    pars = c(lac, mu, K, gam, laa)
  )
  expect_silent(rates <- DAISIE:::update_rates(
    timeval = 0,
    totaltime = 1,
    gam = gam,
    laa = laa,
    lac = lac,
    mu = mu,
    hyper_pars = default_pars$hyper_pars,
    area_pars = default_pars$area_pars,
    dist_pars = default_pars$dist_pars,
    ext_pars = default_pars$ext_pars,
    island_ontogeny = translate_island_ontogeny("beta"),
    sea_level = translate_sea_level("const"),
    extcutoff = 1000.0,
    K = K,
    num_spec = 0,
    num_immigrants = 0,
    mainland_n = 1))
  expect_true(are_rates(rates))
  expected_rates <- list(immig_rate = 0,
                         ext_rate = 0,
                         ana_rate = 0,
                         clado_rate = 0)
  expect_equal(rates, expected_rates)
})

#test_that("update area-dependent rates with sea-level is silent and gives
#correct output"
