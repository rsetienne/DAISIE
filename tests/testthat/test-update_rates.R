context("update_rates")

test_that("update_rates constant rates is silent and gives correct output", {
  timeval <- 0
  totaltime <- 1
  gam <- 0.009
  mu <- 2.0
  laa <- 1.0
  lac <- 2.5
  ddmodel_sim <- 11
  Apars <- NULL
  Epars <- NULL
  island_ontogeny <- translate_island_ontogeny("const")
  extcutoff <- 1000.0
  K <- 3
  island_spec <- c()
  mainland_n <- 1
  t_hor <- 0.5
  set.seed(42)
  expect_silent(rates <- update_rates(
    timeval = timeval,
    totaltime = totaltime,
    gam = gam,
    mu = mu,
    laa = laa,
    lac = lac,
    ddmodel_sim = ddmodel_sim,
    Apars = Apars,
    Epars = Epars,
    island_ontogeny = island_ontogeny,
    extcutoff = extcutoff,
    K = K,
    island_spec = island_spec,
    mainland_n = mainland_n,
    t_hor = t_hor))
  expect_true(are_rates(rates))
  expected_rates <- list(immig_rate = 0.009,
                         ext_rate = 0,
                         ana_rate = 0,
                         clado_rate = 0,
                         ext_rate_max = 0,
                         immig_rate_max = 0.009,
                         clado_rate_max = 0)
  expect_equal(rates, expected_rates)
})


test_that("update area-dependent rates is silent and gives correct output", {
  set.seed(42)
  expect_silent(rates <- update_rates(
    timeval = 0,
    totaltime = 1,
    gam = 0.009,
    mu = 2.0,
    laa = 1.0,
    lac = 2.5,
    ddmodel_sim = 11,
    Apars = create_area_pars(
      max_area = 1.0,
      proportional_peak_t = 0.5,
      peak_sharpness = 1.0,
      total_island_age = 1.0),
    Epars = c(0.5, 10.0),
    island_ontogeny = translate_island_ontogeny("beta"),
    extcutoff = 1000.0,
    K = 3,
    island_spec = c(),
    mainland_n = 1,
    t_hor = 0.5))
  expect_true(are_rates(rates))
  expected_rates <- list(immig_rate = 0,
                         ext_rate = 0,
                         ana_rate = 0,
                         clado_rate = 0,
                         ext_rate_max = 0,
                         immig_rate_max = 0.009,
                         clado_rate_max = 0)
  expect_equal(rates, expected_rates)
})
