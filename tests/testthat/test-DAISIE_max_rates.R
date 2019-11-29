test_that("Yutes at work", {

  totaltime <- 1
  area_pars <- DAISIE::create_area_pars(
    max_area = 5000,
    proportional_peak_t = 0.5,
    peak_sharpness = 1,
    total_island_age = 15,
    sea_level_amplitude = 0,
    sea_level_frequency = 0
  )

  testthat::expect_silent(
    global_max_area_time <- DAISIE:::get_global_max_area_time(
      totaltime = totaltime,
      area_pars = area_pars,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level
    )
  )

})

