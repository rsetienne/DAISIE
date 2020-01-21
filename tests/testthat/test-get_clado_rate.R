context("get_clado_rate")

test_that("use area constant diversity-independent without hyper_pars", {
  default_pars <- create_default_pars(
    island_ontogeny = 0,
    sea_level = 0,
    totaltime = 5,
    area_pars = NULL,
    hyper_pars = NULL,
    dist_pars = NULL,
    ext_pars = NULL,
    pars = c(2, 0, 0, 0, 0)
  )
  expect_silent(
    get_clado_rate(timeval = 0,
                   lac = 2,
                   hyper_pars = default_pars$hyper_pars,
                   area_pars = default_pars$area_pars,
                   dist_pars = default_pars$dist_pars,
                   island_ontogeny = 0,
                   sea_level = 0,
                   num_spec = 0,
                   K = 10
    )
  )
})

test_that("use area constant diversity-independent with hyper_pars", {
  ps_clado_rate <- 0.2
  carr_cap <- 9
  n_species <- 4
  default_pars <- create_default_pars(
    island_ontogeny = 0,
    sea_level = 0,
    totaltime = 5,
    area_pars = NULL,
    hyper_pars = create_hyper_pars(d_0 = 4, x = 3, alpha = 2, beta = 1),
    dist_pars = NULL,
    ext_pars = NULL,
    pars = c(ps_clado_rate, 0, 0, 0, 0)
  )
  created <- get_clado_rate(
    timeval = 5,
    lac = ps_clado_rate,
    hyper_pars = default_pars$hyper_pars,
    area_pars = default_pars$area_pars,
    dist_pars = default_pars$dist_pars,
    island_ontogeny = 0,
    sea_level = 0,
    num_spec = n_species,
    K = carr_cap
  )
  expected <- DAISIE_calc_clade_clado_rate(
    ps_clado_rate = ps_clado_rate,
    n_species = n_species,
    carr_cap = carr_cap
  )
  expect_equal(created, expected)
})

test_that("use area constant diversity-dependent without hyper_pars", {
  ps_clado_rate <- 0.2
  carr_cap <- 9
  n_species <- 4
  default_pars <- create_default_pars(
    island_ontogeny = 0,
    sea_level = 0,
    totaltime = 5,
    area_pars = NULL,
    hyper_pars = NULL,
    dist_pars = NULL,
    ext_pars = NULL,
    pars = c(ps_clado_rate, 0, 0, 0, 0)
  )

  created <- get_clado_rate(
    timeval = 5,
    lac = ps_clado_rate,
    hyper_pars = default_pars$hyper_pars,
    area_pars = default_pars$area_pars,
    dist_pars = default_pars$dist_pars,
    island_ontogeny = 0,
    sea_level = 0,
    num_spec = n_species,
    K = carr_cap
  )
  expected <- DAISIE_calc_clade_clado_rate(
    ps_clado_rate = ps_clado_rate,
    n_species = n_species,
    carr_cap = carr_cap
  )
  expect_equal(created, expected)
})

test_that("use area constant diversity-dependent with hyper_pars", {
  ps_clado_rate <- 0.2
  carr_cap <- 9
  n_species <- 4
  ps_clado_rate <- 0.2
  carr_cap <- 9
  n_species <- 4
  default_pars <- create_default_pars(
    island_ontogeny = 0,
    sea_level = 0,
    totaltime = 5,
    area_pars = NULL,
    hyper_pars = NULL,
    dist_pars = NULL,
    ext_pars = NULL,
    pars = c(ps_clado_rate, 0, 0, 0, 0)
  )

  created <- get_clado_rate(
    timeval = 5,
    lac = ps_clado_rate,
    hyper_pars = default_pars$hyper_pars,
    area_pars = default_pars$area_pars,
    dist_pars = default_pars$dist_pars,
    island_ontogeny = 0,
    sea_level = 0,
    num_spec = n_species,
    K = carr_cap
  )
  expected <- DAISIE_calc_clade_clado_rate(
    ps_clado_rate = ps_clado_rate,
    n_species = n_species,
    carr_cap = carr_cap
  )
  expect_equal(created, expected)
})

test_that("use area constant diversity-independent with hyper_pars", {
  ps_clado_rate <- 0.2
  carr_cap <- Inf
  n_species <- 4
  default_pars <- create_default_pars(
    island_ontogeny = 0,
    sea_level = 0,
    totaltime = 5,
    area_pars = NULL,
    hyper_pars = NULL,
    dist_pars = NULL,
    ext_pars = NULL,
    pars = c(ps_clado_rate, 0, 0, 0, 0)
  )
  created <- get_clado_rate(
    timeval = 5,
    lac = ps_clado_rate,
    hyper_pars = default_pars$hyper_pars,
    area_pars = default_pars$area_pars,
    dist_pars = default_pars$dist_pars,
    island_ontogeny = 0,
    sea_level = 0,
    num_spec = n_species,
    K = carr_cap
  )
  expected <- DAISIE_calc_clade_clado_rate(
    ps_clado_rate = ps_clado_rate,
    n_species = n_species,
    carr_cap = Inf
  )
  expect_equal(created, expected)
})

test_that("use area variable (ontogeny) diversity-dependent without
          hyper_pars",{
            ps_clado_rate <- 0.2
            carr_cap <- 9
            n_species <- 4
            default_pars <- create_default_pars(
              island_ontogeny = 1,
              sea_level = 0,
              totaltime = 5,
              area_pars = create_area_pars(
                max_area = 1,
                proportional_peak_t = 0.5,
                peak_sharpness = 1,
                total_island_age = 10,
                sea_level_amplitude = 0,
                sea_level_frequency = 0
              ),
              hyper_pars = NULL,
              dist_pars = NULL,
              ext_pars = NULL,
              pars = c(ps_clado_rate, 0, 0, 0, 0)
            )
            created <- get_clado_rate(
              timeval = 5,
              lac = ps_clado_rate,
              hyper_pars = default_pars$hyper_pars,
              area_pars = default_pars$area_pars,
              dist_pars = default_pars$dist_pars,
              island_ontogeny = 1,
              sea_level = 0,
              num_spec = n_species,
              K = carr_cap
            )
            expected <- DAISIE_calc_clade_clado_rate(
              ps_clado_rate = ps_clado_rate,
              n_species = n_species,
              carr_cap = carr_cap
            )
            expect_equal(created, expected)
})

test_that("use area variable (ontogeny) diversity-dependent with
          hyper_pars",{
            ps_clado_rate <- 0.2
            default_pars <- create_default_pars(
              island_ontogeny = 0,
              sea_level = 0,
              totaltime = 5,
              area_pars = create_area_pars(1, 0.5, 1, 10, 0, 0),
              hyper_pars = create_hyper_pars(
                d_0 = 4,
                x = 3,
                alpha = 2,
                beta = 1
              ),
              dist_pars = NULL,
              ext_pars = NULL,
              pars = c(ps_clado_rate, 0, 0, 0, 0)
            )
            ps_clado_rate <- 0.2
            carr_cap <- 9
            n_species <- 4
            created <- get_clado_rate(
              timeval = 5,
              lac = ps_clado_rate,
              hyper_pars = default_pars$hyper_pars,
              area_pars = default_pars$area_pars,
              dist_pars = default_pars$dist_pars,
              island_ontogeny = 1,
              sea_level = 0,
              num_spec = n_species,
              K = carr_cap
            )
            expected <- DAISIE_calc_clade_clado_rate(
              ps_clado_rate = ps_clado_rate,
              n_species = n_species,
              carr_cap = carr_cap
            )
            expect_equal(created, expected)
          })

test_that("use area variable (sea-level) diversity-dependent without
          hyper_pars", {
            ps_clado_rate <- 0.2
            default_pars <- create_default_pars(
              island_ontogeny = 0,
              sea_level = 1,
              totaltime = 5,
              area_pars = create_area_pars(1, 0, 0, 10, 1, 1),
              hyper_pars = create_hyper_pars(
                d_0 = 0,
                x = 0,
                alpha = 0,
                beta = 0
              ),
              dist_pars = NULL,
              ext_pars = NULL,
              pars = c(ps_clado_rate, 0, 0, 0, 0)
            )
            ps_clado_rate <- 0.2
            carr_cap <- 9
            n_species <- 4
            created <- get_clado_rate(
              timeval = 5,
              lac = ps_clado_rate,
              hyper_pars = default_pars$hyper_pars,
              area_pars = default_pars$area_pars,
              dist_pars = default_pars$dist_pars,
              island_ontogeny = 0,
              sea_level = 1,
              num_spec = n_species,
              K = carr_cap
            )
            expected <- DAISIE_calc_clade_clado_rate(
              ps_clado_rate = ps_clado_rate,
              n_species = n_species,
              carr_cap = carr_cap
            )
            expect_equal(created, expected)
          })

test_that("use area variable (sea-level) diversity-dependent with
          hyper_pars", {
            ps_clado_rate <- 0.2
            default_pars <- create_default_pars(
              island_ontogeny = 0,
              sea_level = 1,
              totaltime = 5,
              area_pars = create_area_pars(1, 0, 0, 10, 1, 1),
              hyper_pars = create_hyper_pars(
                d_0 = 1,
                x = 3,
                alpha = 2,
                beta = 1
              ),
              dist_pars = NULL,
              ext_pars = NULL,
              pars = c(ps_clado_rate, 0, 0, 0, 0)
            )
            ps_clado_rate <- 0.2
            carr_cap <- 9
            n_species <- 4
            created <- get_clado_rate(
              timeval = 5,
              lac = ps_clado_rate,
              hyper_pars = default_pars$hyper_pars,
              area_pars = default_pars$area_pars,
              dist_pars = default_pars$dist_pars,
              island_ontogeny = 0,
              sea_level = 1,
              num_spec = n_species,
              K = carr_cap
            )
            expected <- DAISIE_calc_clade_clado_rate(
              ps_clado_rate = ps_clado_rate,
              n_species = n_species,
              carr_cap = carr_cap
            )
            expect_equal(created, expected)
          })

test_that("use area variable (ontogeny and sea-level) diversity-dependent
          without hyper_pars", {
            ps_clado_rate <- 0.2
            default_pars <- create_default_pars(
              island_ontogeny = 0,
              sea_level = 1,
              totaltime = 5,
              area_pars = create_area_pars(1, 0.5, 1, 10, 1, 1),
              hyper_pars = create_hyper_pars(
                d_0 = 0,
                x = 0,
                alpha = 0,
                beta = 0
              ),
              dist_pars = NULL,
              ext_pars = NULL,
              pars = c(ps_clado_rate, 0, 0, 0, 0)
            )
            ps_clado_rate <- 0.2
            carr_cap <- 9
            n_species <- 4
            created <- get_clado_rate(
              timeval = 5,
              lac = ps_clado_rate,
              hyper_pars = default_pars$hyper_pars,
              area_pars = default_pars$area_pars,
              dist_pars = default_pars$dist_pars,
              island_ontogeny = 1,
              sea_level = 1,
              num_spec = n_species,
              K = carr_cap
            )
            expected <- DAISIE_calc_clade_clado_rate(
              ps_clado_rate = ps_clado_rate,
              n_species = n_species,
              carr_cap = carr_cap
            )
            expect_equal(created, expected)
          })

test_that("use area variable (ontogeny and sea-level) diversity-dependent
          with hyper_pars", {
            ps_clado_rate <- 0.2
            default_pars <- create_default_pars(
              island_ontogeny = 0,
              sea_level = 1,
              totaltime = 5,
              area_pars = create_area_pars(1, 0.5, 1, 10, 1, 1),
              hyper_pars = create_hyper_pars(
                d_0 = 1,
                x = 3,
                alpha = 2,
                beta = 1
              ),
              dist_pars = NULL,
              ext_pars = NULL,
              pars = c(ps_clado_rate, 0, 0, 0, 0)
            )
            ps_clado_rate <- 0.2
            carr_cap <- 9
            n_species <- 4
            created <- get_clado_rate(
              timeval = 5,
              lac = ps_clado_rate,
              hyper_pars = default_pars$hyper_pars,
              area_pars = default_pars$area_pars,
              dist_pars = default_pars$dist_pars,
              island_ontogeny = 1,
              sea_level = 1,
              num_spec = n_species,
              K = carr_cap
            )
            expected <- DAISIE_calc_clade_clado_rate(
              ps_clado_rate = ps_clado_rate,
              n_species = n_species,
              carr_cap = carr_cap
            )
            expect_equal(created, expected)
          })
