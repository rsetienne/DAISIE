context("get_dynamic_t_hor")

test_that("use ontogeny timeval < t_hor", {
  t_hor <- 5
  timeval <- 1
  totaltime <- 10
  ext <- 0.5
  area_pars <- DAISIE::create_area_pars(
    max_area = 5000,
    proportional_peak_t = 0.5,
    peak_sharpness = 1,
    total_island_age = 15,
    sea_level_amplitude = 0,
    sea_level_frequency = 0
  )
  island_ontogeny <- 1
  sea_level <- 0
  ext_multiplier <- 1000

  expect_silent(
    dynamic_t_hor <- DAISIE:::get_dynamic_t_hor(
      t_hor = t_hor,
      timeval = timeval,
      totaltime = totaltime,
      ext = ext,
      area_pars = area_pars,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      ext_multiplier = ext_multiplier
    )
  )

})
test_that("use t_hor NULL", {
  t_hor <- NULL
  timeval <- 1
  totaltime <- 10
  ext <- 0.5
  area_pars <- DAISIE::create_area_pars(
    max_area = 5000,
    proportional_peak_t = 0.5,
    peak_sharpness = 1,
    total_island_age = 15,
    sea_level_amplitude = 0,
    sea_level_frequency = 0
  )
  island_ontogeny <- 1
  sea_level <- 0
  ext_multiplier <- 1000

  expect_silent(
    dynamic_t_hor <- DAISIE:::get_dynamic_t_hor(
      t_hor = t_hor,
      timeval = timeval,
      totaltime = totaltime,
      ext = ext,
      area_pars = area_pars,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      ext_multiplier = ext_multiplier
    )
  )

})

test_that("use ontogeny timeval > t_hor", {
  t_hor <- 5
  timeval <- 6
  totaltime <- 10
  ext <- 0.5
  area_pars <- DAISIE::create_area_pars(
    max_area = 5000,
    proportional_peak_t = 0.5,
    peak_sharpness = 1,
    total_island_age = 15,
    sea_level_amplitude = 0,
    sea_level_frequency = 0
  )
  island_ontogeny <- 1
  sea_level <- 0
  ext_multiplier <- 1000

  expect_silent(
    dynamic_t_hor <- DAISIE:::get_dynamic_t_hor(
      t_hor = t_hor,
      timeval = timeval,
      totaltime = totaltime,
      ext = ext,
      area_pars = area_pars,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      ext_multiplier = ext_multiplier
    )
  )

})

test_that("abuse", {
  t_hor <- 5
  timeval <- 1
  totaltime <- 10
  ext <- 0.5
  area_pars <- DAISIE::create_area_pars(
    max_area = 5000,
    proportional_peak_t = 0.5,
    peak_sharpness = 1,
    total_island_age = 15,
    sea_level_amplitude = 0,
    sea_level_frequency = 0
  )
  island_ontogeny <- 1
  sea_level <- 0
  ext_multiplier <- 1000
  expect_error(
    dynamic_t_hor <- DAISIE:::get_dynamic_t_hor(
      t_hor = "nonsense",
      timeval = timeval,
      totaltime = totaltime,
      ext = ext,
      area_pars = area_pars,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      ext_multiplier = ext_multiplier
    )
  )
})
