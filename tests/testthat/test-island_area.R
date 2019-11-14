context("island_area")

test_that("is valid island area", {
  area <- island_area(timeval = 2,
                      area_pars = create_area_pars(
                      max_area = 10,
                      proportional_peak_t = 0.5,
                      peak_sharpness = 1,
                      total_island_age = 5,
                      sea_level_amplitude = 5,
                      sea_level_frequency = 10),
                      island_ontogeny = translate_island_ontogeny("beta"),
                      sea_level = translate_sea_level("const")
                      )
  expect_true(is.numeric(area) && area >= 0)
})
