context("island_area")

test_that("is valid island area with constant area", {
  area <- island_area(timeval = 2,
                      area_pars <- create_area_pars(
                        max_area = 1,
                        current_area = 1,
                        proportional_peak_t = 0,
                        total_island_age = 5,
                        sea_level_amplitude = 0,
                        sea_level_frequency = 0,
                        island_gradient_angle = 0),
                      island_ontogeny = translate_island_ontogeny("const"),
                      sea_level = translate_sea_level("const")
                      )
  expect_true(is.numeric(area) && area >= 0)
})

test_that("is valid island area with ontogeny", {
  area_pars <- create_area_pars(
    max_area = 10,
    current_area = 5,
    proportional_peak_t = 0.5,
    total_island_age = 5,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0
  )
  peak <- DAISIE:::calc_peak(totaltime = 4, area_pars = area_pars)

  area <- island_area(timeval = 2,
                      area_pars = area_pars,
                      peak = peak,
                      island_ontogeny = translate_island_ontogeny("beta"),
                      sea_level = translate_sea_level("const")
  )
  expect_true(is.numeric(area) && area >= 0)
})

test_that("is valid island area with sea level", {
  area <- island_area(timeval = 2,
                      area_pars = create_area_pars(
                        max_area = 1000,
                        current_area = 1,
                        proportional_peak_t = 0,
                        total_island_age = 10,
                        sea_level_amplitude = 60,
                        sea_level_frequency = 3,
                        island_gradient_angle = 85),
                      island_ontogeny = translate_island_ontogeny("const"),
                      sea_level = translate_sea_level("sine")
  )
  expect_true(is.numeric(area) && area >= 0)
})

test_that("is valid island area with ontogeny and sea level", {
  area_pars <- create_area_pars(
    max_area = 10,
    current_area = 2,
    proportional_peak_t = 0.5,
    total_island_age = 5,
    sea_level_amplitude = 2,
    sea_level_frequency = 10,
    island_gradient_angle = 85
  )
  peak <- DAISIE:::calc_peak(totaltime = 4, area_pars = area_pars)
  area <- island_area(timeval = 2,
                      area_pars = area_pars,
                      peak = peak,
                      island_ontogeny = translate_island_ontogeny("beta"),
                      sea_level = translate_sea_level("sine")
  )
  expect_true(is.numeric(area) && area >= 0)
})

test_that("abuse island area with constant area", {
  expect_warning(island_area(timeval = 2,
                             area_pars = create_area_pars(
                               max_area = 10,
                               current_area = 1,
                               proportional_peak_t = 0.5,
                               total_island_age = 5,
                               sea_level_amplitude = 5,
                               sea_level_frequency = 10,
                               island_gradient_angle = 45),
                             island_ontogeny = translate_island_ontogeny("const"),
                             sea_level = translate_sea_level("const")),
                 "Constant island area requires a maximum area of 1.")
})
