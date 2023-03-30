test_that("use ontogeny", {

  total_time <- 1
  area_pars <- create_area_pars(
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
  peak <- calc_peak(total_time = total_time, area_pars = area_pars)
  testthat::expect_silent(
    global_max_area <- get_global_max_area(
      total_time = total_time,
      area_pars = area_pars,
      peak = peak,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level
    )
  )
  total_time <- 15
  # Gets actual peak for entire curve
  testthat::expect_equal(
    global_max_area <- get_global_max_area(
      total_time = total_time,
      area_pars = area_pars,
      peak = peak,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level
    ),
    expected = area_pars$max_area
  )
})

test_that("use sea level", {

  total_time <- 1
  area_pars <- create_area_pars(
    max_area = 1000,
    current_area = 1,
    proportional_peak_t = 0,
    total_island_age = 15,
    sea_level_amplitude = 50,
    sea_level_frequency = 10,
    island_gradient_angle = 45
  )

  island_ontogeny <- 0
  sea_level <- 1
  testthat::expect_silent(
    global_max_area <- get_global_max_area(
      total_time = total_time,
      area_pars = area_pars,
      peak = peak,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level
    )
  )
})
