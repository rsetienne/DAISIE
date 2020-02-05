context("get_global_max_area_time")

test_that("use ontogeny", {

  totaltime <- 1
  area_pars <- DAISIE::create_area_pars(
    max_area = 5000,
    proportional_peak_t = 0.5,
    peak_sharpness = 1,
    total_island_age = 15,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0
  )
  island_ontogeny <- 1
  sea_level <- 0
  testthat::expect_silent(
    global_max_area_time <- DAISIE:::get_global_max_area_time(
      totaltime = totaltime,
      area_pars = area_pars,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level
    )
  )

})

test_that("use sea level", {

  totaltime <- 1
  area_pars <- DAISIE::create_area_pars(
    max_area = 1000,
    proportional_peak_t = 0,
    peak_sharpness = 0,
    total_island_age = 15,
    sea_level_amplitude = 50,
    sea_level_frequency = 10,
    island_gradient_angle = 45
  )

  island_ontogeny <- 0
  sea_level <- 1
  testthat::expect_silent(
    global_max_area_time <- DAISIE:::get_global_max_area_time(
      totaltime = totaltime,
      area_pars = area_pars,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level
    )
  )
})
