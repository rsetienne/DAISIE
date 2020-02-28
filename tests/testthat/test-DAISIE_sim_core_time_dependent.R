context("DAISIE_sim_core_time_dependent")

test_that("Ontogeny oceanic should run silent IW", {
  set.seed(234567890)
  expect_silent(
    DAISIE:::DAISIE_sim_core_time_dependent(
      time = 10,
      mainland_n = 100,
      pars = c(0.0001, 2.2, 0.005, 0.001, 1),
      area_pars = create_area_pars(
        max_area = 5000,
        proportional_peak_t = 0.5,
        peak_sharpness = 1,
        total_island_age = 15,
        sea_level_amplitude = 0,
        sea_level_frequency = 0,
        island_gradient_angle = 0
      ),
      ext_pars = c(1, 100),
      island_ontogeny = "beta",
      sea_level = "const",
      extcutoff = 20 #1000
    )
  )
})

test_that("Ontogeny oceanic should run silent CS", {
  set.seed(420)
  expect_silent(
    DAISIE:::DAISIE_sim_core_time_dependent(
      time = 10,
      mainland_n = 1,
      pars = c(0.0001, 2.2, 0.005, 0.001, 1),
      area_pars = create_area_pars(
        max_area = 5000,
        proportional_peak_t = 0.5,
        peak_sharpness = 1,
        total_island_age = 15,
        sea_level_amplitude = 0,
        sea_level_frequency = 0,
        island_gradient_angle = 0
      ),
      ext_pars = c(1, 100),
      island_ontogeny = "beta",
      sea_level = "const"
    )
  )
})

test_that("Ontogeny oceanic with sea level should run silent CS", {
  set.seed(439)
  expect_silent(
    out <- DAISIE:::DAISIE_sim_core_time_dependent(
      time = 10,
      mainland_n = 1,
      pars = c(0.00001, 2.2, 0.005, 0.06, 1),
      area_pars = create_area_pars(
        max_area = 5000,
        proportional_peak_t = 0.5,
        peak_sharpness = 1,
        total_island_age = 15,
        sea_level_amplitude = 60,
        sea_level_frequency = 10,
        island_gradient_angle = 30
      ),
      ext_pars = c(1, 100),
      island_ontogeny = "beta",
      sea_level = "const"
    )
  )
  expect_equal(out$branching_times, c(10, 0.17840243993784999))
})

test_that("all species extinct if island dead", {
  ontogeny_sim <- DAISIE:::DAISIE_sim_core_time_dependent(
    time = 10,
    mainland_n = 1000,
    pars = c(0.0001, 2.2, 0.005, 0.001, 1),
    area_pars = create_area_pars(
      max_area = 5000,
      proportional_peak_t = 0.5,
      peak_sharpness = 1,
      total_island_age = 10,
      sea_level_amplitude = 0,
      sea_level_frequency = 0,
      island_gradient_angle = 0
    ),
    ext_pars = c(1, 100),
    island_ontogeny = "beta",
    sea_level = "const",
    extcutoff = 20 #1000
  )
  last_entry <- ontogeny_sim$stt_table[nrow(ontogeny_sim$stt_table), ]
  expect_true(last_entry[1] == 0)
  expect_true(last_entry[2] == 0)
  expect_true(last_entry[3] == 0)
  expect_true(last_entry[4] == 0)
})

test_that("!is.null(area_pars) && island_ontogeny == 'const'", {
  expect_error(DAISIE:::DAISIE_sim_core_time_dependent(time = 1,
                                                       mainland_n = 100,
                                                       pars = c(2, 2, 20, 0, 1),
                                                       island_ontogeny = 0,
                                                       sea_level = 0,
                                                       area_pars = create_area_pars(
                                                         max_area = 1,
                                                         proportional_peak_t = 1,
                                                         peak_sharpness = 1,
                                                         total_island_age = 1,
                                                         sea_level_amplitude = 0,
                                                         sea_level_frequency = 0,
                                                         island_gradient_angle = 0
                                                       )
  ), regexp = "area_pars specified for constant island_ontogeny and sea_level.
         Run DAISIE_sim_constant_rate instead.")
})

test_that("(is.null(ext_pars) || is.null(area_pars)) &&
          (island_ontogeny != 0 || sea_level != 0)", {

            time <- 10
            mainland_n <- 1000
            pars <- c(0.0001, 2.2, 0.005, 0.001, 1)
            area_pars <- NULL
            ext_pars <- NULL
            island_ontogeny <- 1
            sea_level <- 1

            expect_error(
              DAISIE:::DAISIE_sim_core_time_dependent(
                time = time,
                mainland_n = mainland_n,
                pars = pars,
                area_pars = area_pars,
                ext_pars = ext_pars,
                island_ontogeny = island_ontogeny,
                sea_level = sea_level
              ), regexp =
                "Island ontogeny and/or sea level specified but area parameters
    and/or extinction parameters not available. Please either set
    island_ontogeny and sea_level to NULL, or specify area_pars and ext_pars."
            )
          })


