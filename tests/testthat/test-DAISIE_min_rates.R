context("get_global_min_area")

test_that("use ontogeny", {

  totaltime <- 1
  area_pars <- DAISIE::create_area_pars(
    max_area = 5000,
    current_area = 2500,
    proportional_peak_t = 0.5,
    total_island_age = 15,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0
  )
  island_ontogeny <- 1
  sea_level <- 0
  peak <- DAISIE:::calc_peak(totaltime = totaltime, area_pars = area_pars)
  testthat::expect_silent(
    global_min_area <- DAISIE:::get_global_min_area(
      totaltime = totaltime,
      area_pars = area_pars,
      peak = peak,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level
    )
  )
})

test_that("use sea level", {

  totaltime <- 1
  area_pars <- DAISIE::create_area_pars(
    max_area = 1000,
    current_area = 1,
    proportional_peak_t = 0,
    total_island_age = 15,
    sea_level_amplitude = 50,
    sea_level_frequency = 10,
    island_gradient_angle = 85
  )

  island_ontogeny <- 0
  sea_level <- 1
  testthat::expect_silent(
    global_min_area <- DAISIE:::get_global_min_area(
      totaltime = totaltime,
      area_pars = area_pars,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level
    )
  )
  testthat::expect_equal(global_min_area, 569.74327591692145)
})
