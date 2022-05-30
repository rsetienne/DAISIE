test_that("update_rates constant rates is silent and gives correct output", {
  timeval <- 0
  total_time <- 1
  gam <- 0.009
  mu <- 2.0
  laa <- 1.0
  lac <- 2.5
  K <- 3
  area_pars <- create_area_pars(max_area = 1,
                                current_area = 1,
                                proportional_peak_t = 0,
                                total_island_age = 0,
                                sea_level_amplitude = 0,
                                sea_level_frequency = 0,
                                island_gradient_angle = 0)
  hyper_pars <- create_hyper_pars(d = 0, x = 0)
  island_ontogeny <- translate_island_ontogeny("const")
  sea_level <- translate_sea_level("const")
  extcutoff <- 1000.0
  num_spec <- 0
  num_immigrants <- 0
  mainland_n <- 1
  set.seed(42)
  expect_silent(rates <- update_rates(
    timeval = timeval,
    total_time = total_time,
    gam = gam,
    laa = laa,
    lac = lac,
    mu = mu,
    hyper_pars = hyper_pars,
    area_pars = area_pars,
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
  area_pars <- create_area_pars(
    max_area = 1.0,
    current_area = 0.5,
    proportional_peak_t = 0.5,
    total_island_age = 1.2,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0)
  hyper_pars <- create_hyper_pars(d = 0.2, x = 0.1)
  peak <- calc_peak(total_time = 1, area_pars = area_pars)
  expect_silent(rates <- DAISIE:::update_rates(
    timeval = 0,
    total_time = 1,
    gam = gam,
    laa = laa,
    lac = lac,
    mu = mu,
    hyper_pars = hyper_pars,
    area_pars = area_pars,
    island_ontogeny = translate_island_ontogeny("beta"),
    sea_level = translate_sea_level("const"),
    extcutoff = 1000.0,
    K = K,
    num_spec = 0,
    num_immigrants = 0,
    mainland_n = 1,
    peak = peak))
  expect_true(are_rates(rates))
  expected_rates <- list(immig_rate = 0,
                         ext_rate = 0,
                         ana_rate = 0,
                         clado_rate = 0)
  expect_equal(rates, expected_rates)
})

#test_that("update area-dependent rates with sea-level is silent and gives
#correct output"
