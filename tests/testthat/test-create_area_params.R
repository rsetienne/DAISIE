context("create_area_params")

test_that("minimal use", {
  expect_silent(
    create_area_params(
      max_area = 10,
      proportional_peak_t = 0.5,
      peak_sharpness = 1,
      total_island_age = 5
    )
  )

})
