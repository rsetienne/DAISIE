context("are_area_pars")

test_that("minimal use", {
  expect_true(
    are_area_pars(
      create_area_pars(
        max_area = 10,
        current_area = 1,
        proportional_peak_t = 0.5,
        total_island_age = 5,
        sea_level_amplitude = 5,
        sea_level_frequency = 10,
        island_gradient_angle = 45)))
})
