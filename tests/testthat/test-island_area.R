context("island_area")

test_that("is valid island area", {
  area <- island_area(timeval = 2,
                      Apars = create_area_params(
                      max_area = 10,
                      proportional_peak_t = 0.5,
                      peak_sharpness = 1,
                      total_island_age = 5),island_ontogeny = "quadratic"
                      )
  expect_true(is.numeric(area) && area >= 0)
})
