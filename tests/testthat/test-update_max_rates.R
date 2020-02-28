context("update_max_rates")

test_that("update_max_rates constant rates is silent and gives correct output", {
  timeval <- 0
  totaltime <- 1
  gam <- 0.009
  laa <- 1.0
  lac <- 2.5
  mu <- 2.5
  default_pars <-
    create_default_pars(
      island_ontogeny = 0,
      sea_level = 0,
      area_pars = create_area_pars(
        max_area = 1,
        proportional_peak_t = 0,
        peak_sharpness = 0,
        total_island_age = 1,
        sea_level_amplitude = 0,
        sea_level_frequency = 0,
        island_gradient_angle = 0
      ),
      hyper_pars = create_hyper_pars(0, 0, 0, 0),
      dist_pars = create_dist_pars(1),
      ext_pars = c(2),
      totaltime = totaltime,
      pars = c(0, 0, 0, 0, 0)
    )
  island_ontogeny <- translate_island_ontogeny("const")
  sea_level <- translate_sea_level("const")
  extcutoff <- 1000.0
  K <- 3
  num_spec <- 0
  num_immigrants <- 0
  mainland_n <- 1
  global_min_area_time <- get_global_min_area_time(totaltime = totaltime,
                                                   area_pars = default_pars$area_pars,
                                                   island_ontogeny = island_ontogeny,
                                                   sea_level = sea_level)
  global_max_area_time <- get_global_max_area_time(totaltime = totaltime,
                                                   area_pars = default_pars$area_pars,
                                                   island_ontogeny = island_ontogeny,
                                                   sea_level = sea_level)
  set.seed(42)
  expect_silent(rates <- update_max_rates(
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
    mainland_n = mainland_n,
    global_min_area_time = global_min_area_time,
    global_max_area_time = global_max_area_time))
  expect_true(are_max_rates(rates))
  expected_rates <- list(
    ext_max_rate = 0,
    immig_max_rate = 0.009,
    ana_max_rate = 0,
    clado_max_rate = 0
  )
  expect_equal(rates, expected_rates)
})


test_that("update area-dependent max rates is silent and gives correct output", {
  set.seed(42)

  default_pars <-
    create_default_pars(
      island_ontogeny = 1,
      sea_level = 0,
      area_pars = create_area_pars(
        max_area = 1.0,
        proportional_peak_t = 0.5,
        peak_sharpness = 1.0,
        total_island_age = 1.0,
        sea_level_amplitude = 0,
        sea_level_frequency = 0,
        island_gradient_angle = 0
      ),
      hyper_pars = create_hyper_pars(0, 0, 0, 0),
      dist_pars = create_dist_pars(1),
      ext_pars = c(0.5, 10.0),
      totaltime = totaltime,
      pars = c(0, 0, 0, 0, 0)
    )
  global_min_area_time <- get_global_min_area_time(totaltime = 1,
                                                   area_pars = default_pars$area_pars,
                                                   island_ontogeny = 1,
                                                   sea_level = 0)
  global_max_area_time <- get_global_max_area_time(totaltime = 1,
                                                   area_pars = default_pars$area_pars,
                                                   island_ontogeny = 1,
                                                   sea_level = 0)



  expect_silent(rates <- DAISIE:::update_max_rates(
    timeval = 0,
    totaltime = 1,
    gam = 0.009,
    laa = 1.0,
    lac = 2.5,
    mu = 2.5,
    hyper_pars = default_pars$hyper_pars,
    area_pars = default_pars$area_pars,
    dist_pars = default_pars$dist_pars,
    ext_pars = c(0.5, 10.0),
    island_ontogeny = translate_island_ontogeny("beta"),
    sea_level = translate_sea_level("const"),
    extcutoff = 1000.0,
    K = 3,
    num_spec = 0,
    num_immigrants = 0,
    mainland_n = 1,
    global_min_area_time = global_min_area_time,
    global_max_area_time = global_max_area_time))
  expect_true(are_max_rates(rates))
  expected_rates <- list(
    ext_max_rate = 0,
    immig_max_rate = 0.009,
    ana_max_rate = 0,
    clado_max_rate = 0
  )
  expect_equal(rates, expected_rates)
})

#test_that("update area-dependent rates with sea-level is silent and gives
#correct output"
