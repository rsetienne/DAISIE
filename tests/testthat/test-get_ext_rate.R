context("get_ext_rate")

test_that("use area constant diversity-independent without hyper_pars", {
  carr_cap <- 10
  ps_ext_rate <- 2
  n_species <- 5
  n_mainland_species <- 2
  default_pars <- create_default_pars(
    island_ontogeny = 1,
    sea_level = 0,
    totaltime = 5,
    area_pars = create_area_pars(
      max_area = 1,
      proportional_peak_t = 0,
      peak_sharpness = 0,
      total_island_age = 0,
      sea_level_amplitude = 0,
      sea_level_frequency = 0
    ),
    hyper_pars = NULL,
    dist_pars = NULL,
    ext_pars = c(),
    pars = c(0, ps_ext_rate, 0, 0, 0)
  )

  expect_silent(ext_rate <- DAISIE:::get_ext_rate(
    timeval = 0,
    ext_pars = default_pars$ext_pars,
    hyper_pars = default_pars$hyper_pars,
    area_pars = default_pars$area_pars,
    island_ontogeny = translate_island_ontogeny("const"),
    sea_level = translate_sea_level("const"),
    extcutoff = 100,
    num_spec = 0,
    K = 10)
  )
  expect_true(is.numeric(ext_rate))
  expected <- DAISIE_calc_clade_ext_rate(
    ps_ext_rate = ps_ext_rate,
    n_species = n_species
  )
  created <- get_ext_rate(
    timeval = 10,
    ext_pars = default_pars$ext_pars,
    hyper_pars = default_pars$hyper_pars,
    area_pars =  default_pars$area_pars,
    island_ontogeny = translate_island_ontogeny("const"),
    sea_level = translate_sea_level("const"),
    extcutoff = 100,
    num_spec = n_species,
    K = carr_cap
  )
  expect_equal(expected, created)
})



test_that("use area constant diversity-independent with hyper_pars", {
  ps_ext_rate <- 2
  default_pars <- create_default_pars(
    island_ontogeny = 1,
    sea_level = 0,
    totaltime = 5,
    area_pars = create_area_pars(
      max_area = 1,
      proportional_peak_t = 0,
      peak_sharpness = 0,
      total_island_age = 0,
      sea_level_amplitude = 0,
      sea_level_frequency = 0
    ),
    hyper_pars = create_hyper_pars(d_0 = 4, x = 3, alpha = 2, beta = 1),
    dist_pars = NULL,
    ext_pars = NULL,
    pars = c(0, ps_ext_rate, 0, 0, 0)
  )
  expect_silent(ext_rate <- DAISIE:::get_ext_rate(
    timeval = 0,
    ext_pars = default_pars$ext_pars,
    hyper_pars = default_pars$hyper_pars,
    area_pars = default_pars$area_pars,
    island_ontogeny = translate_island_ontogeny("const"),
    sea_level = translate_sea_level("const"),
    extcutoff = 100,
    num_spec = 0,
    K = 10)
  )
  expect_true(is.numeric(ext_rate))
})


test_that("use area variable (ontogeny) diversity-independent without
          hyper_pars", {
            ps_ext_rate <- 2

            default_pars <- create_default_pars(
              island_ontogeny = 1,
              sea_level = 0,
              totaltime = 5,
              area_pars = create_area_pars(
                max_area = 1000,
                proportional_peak_t = 0.5,
                peak_sharpness = 1,
                total_island_age = 15,
                sea_level_amplitude = 0,
                sea_level_frequency = 0
              ),
              hyper_pars = NULL,
              dist_pars = NULL,
              ext_pars = c(1, 10),
              pars = c(0, ps_ext_rate, 0, 0, 0)
            )

            expect_silent(ext_rate <- DAISIE:::get_ext_rate(
              timeval = 5,
              hyper_pars = NULL,
              area_pars = create_area_pars(1000, 0.5, 1, 15, 0, 0),
              ext_pars = c(1, 10),
              island_ontogeny = translate_island_ontogeny("beta"),
              sea_level = translate_sea_level("const"),
              extcutoff = 1000,
              num_spec = 10,
              K = 20
            )
            )
            expect_true(is.numeric(ext_rate))
          })



test_that("use area variable (sea-level) diversity-independent without
          hyper_pars", {
            ps_ext_rate <- 2

            default_pars <- create_default_pars(
              island_ontogeny = 0,
              sea_level = 1,
              totaltime = 5,
              area_pars = create_area_pars(
                max_area = 10,
                proportional_peak_t = 0,
                peak_sharpness = 0,
                total_island_age = 11,
                sea_level_amplitude = 1,
                sea_level_frequency = 1
              ),
              hyper_pars = NULL,
              dist_pars = NULL,
              ext_pars = c(1, 10),
              pars = c(0, ps_ext_rate, 0, 0, 0)
            )

            expect_silent(ext_rate <- DAISIE:::get_ext_rate(
              timeval = 5,
              hyper_pars = default_pars$hyper_pars,
              area_pars = default_pars$area_pars,
              ext_pars = default_pars$ext_pars,
              island_ontogeny = translate_island_ontogeny("const"),
              sea_level = translate_sea_level("sine"),
              extcutoff = 100,
              num_spec = 10,
              K = 20
            )
            )
          })


test_that("use area variable (ontogeny and sea-level) diversity-independent
          without hyper_pars", {
            ps_ext_rate <- 2

            default_pars <- create_default_pars(
              island_ontogeny = 0,
              sea_level = 1,
              totaltime = 5,
              area_pars = create_area_pars(
                max_area = 1000,
                proportional_peak_t = 0.3,
                peak_sharpness = 1,
                total_island_age = 11,
                sea_level_amplitude = 1,
                sea_level_frequency = 1
              ),
              hyper_pars = NULL,
              dist_pars = NULL,
              ext_pars = c(1, 10),
              pars = c(0, ps_ext_rate, 0, 0, 0)
            )

            expect_silent(ext_rate <- DAISIE:::get_ext_rate(
              timeval = 5,
              hyper_pars = default_pars$hyper_pars,
              area_pars = default_pars$area_pars,
              ext_pars = default_pars$ext_pars,
              island_ontogeny = translate_island_ontogeny("beta"),
              sea_level = translate_sea_level("sine"),
              extcutoff = 100,
              num_spec = 10,
              K = 20
            )
            )
          })

test_that("use area variable (ontogeny and sea-level) diversity-independent
          with hyper_pars", {
            ps_ext_rate <- 2

            default_pars <- create_default_pars(
              island_ontogeny = 0,
              sea_level = 1,
              totaltime = 5,
              area_pars = create_area_pars(
                max_area = 1000,
                proportional_peak_t = 0.3,
                peak_sharpness = 1,
                total_island_age = 11,
                sea_level_amplitude = 1,
                sea_level_frequency = 1
              ),
              hyper_pars = create_hyper_pars(
                d_0 = 4,
                x = 3,
                alpha = 2,
                beta = 1
              ),
              dist_pars = NULL,
              ext_pars = c(1, 10),
              pars = c(0, ps_ext_rate, 0, 0, 0)
            )

            expect_silent(ext_rate <- DAISIE:::get_ext_rate(
              timeval = 5,
              hyper_pars = default_pars$hyper_pars,
              area_pars = default_pars$area_pars,
              ext_pars = default_pars$ext_pars,
              island_ontogeny = translate_island_ontogeny("beta"),
              sea_level = translate_sea_level("sine"),
              extcutoff = 100,
              num_spec = 10,
              K = Inf
            )
            )
          })
