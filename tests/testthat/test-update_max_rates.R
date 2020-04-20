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
        current_area = 0.5,
        proportional_peak_t = 0,
        total_island_age = 1,
        sea_level_amplitude = 0,
        sea_level_frequency = 0,
        island_gradient_angle = 0
      ),
      hyper_pars = create_hyper_pars(0, 0),
      totaltime = totaltime
    )
  island_ontogeny <- translate_island_ontogeny("const")
  sea_level <- translate_sea_level("const")
  extcutoff <- 1000.0
  K <- 3
  num_spec <- 0
  num_immigrants <- 0
  mainland_n <- 1
  peak <- 1
  Amin <- get_global_min_area(totaltime = totaltime,
                              area_pars = default_pars$area_pars,
                              peak = peak,
                              island_ontogeny = island_ontogeny,
                              sea_level = sea_level)
  Amax <- get_global_max_area(totaltime = totaltime,
                              area_pars = default_pars$area_pars,
                              peak = peak,
                              island_ontogeny = island_ontogeny,
                              sea_level = sea_level)
  set.seed(42)
  expect_silent(rates <- update_max_rates(
    gam = gam,
    laa = laa,
    lac = lac,
    mu = mu,
    hyper_pars = default_pars$hyper_pars,
    extcutoff = extcutoff,
    K = K,
    num_spec = num_spec,
    num_immigrants = num_immigrants,
    mainland_n = mainland_n,
    Amin = Amin,
    Amax = Amax))
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
        max_area = 100,
        current_area = 10,
        proportional_peak_t = 0.5,
        total_island_age = 1.0,
        sea_level_amplitude = 0,
        sea_level_frequency = 0,
        island_gradient_angle = 0
      ),
      hyper_pars = create_hyper_pars(0.2, 0.1),
      totaltime = 1
    )
  peak <- DAISIE:::calc_peak(totaltime = 0.7, area_pars = default_pars$area_pars)
  Amin <- get_global_min_area(totaltime = 0.7,
                              area_pars = default_pars$area_pars,
                              peak = peak,
                              island_ontogeny = 1,
                              sea_level = 0)
  Amax <- get_global_max_area(totaltime = 0.7,
                              area_pars = default_pars$area_pars,
                              peak = peak,
                              island_ontogeny = 1,
                              sea_level = 0)



  expect_silent(rates <- DAISIE:::update_max_rates(
    gam = 0.009,
    laa = 1.0,
    lac = 2.5,
    mu = 2.5,
    hyper_pars = default_pars$hyper_pars,
    extcutoff = 1000.0,
    K = 3,
    num_spec = 0,
    num_immigrants = 0,
    mainland_n = 1,
    Amin = Amin,
    Amax = Amax))
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
