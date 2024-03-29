test_that("minimal use", {
  testthat::expect_silent(
    create_area_pars(
      max_area = 10,
      current_area = 5,
      proportional_peak_t = 0.5,
      total_island_age = 5,
      sea_level_amplitude = 5,
      sea_level_frequency = 10,
      island_gradient_angle = 45
    )
  )
  out <- create_area_pars(
    max_area = 10,
    current_area = 5,
    proportional_peak_t = 0.5,
    total_island_age = 5,
    sea_level_amplitude = 5,
    sea_level_frequency = 10,
    island_gradient_angle = 45
  )
  reference <- list(
    max_area = 10,
    current_area = 5,
    proportional_peak_t = 0.5,
    total_island_age = 5,
    sea_level_amplitude = 5,
    sea_level_frequency = 10,
    island_gradient_angle = 45
  )
  testthat::expect_equal(out, reference)
})

test_that("abuse", {
  testthat::expect_error(
    create_area_pars(
      max_area = 0,
      current_area = 0,
      proportional_peak_t = 0.5,
      total_island_age = 5,
      sea_level_amplitude = 5,
      sea_level_frequency = 10,
      island_gradient_angle = 45
    ), regexp = "max_area > 0 is not TRUE"
  )
})
