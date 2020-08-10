context("DAISIE_sim_core_time_dependent")

test_that("Ontogeny oceanic should run silent IW", {
  totaltime <- 10
  mainland_n <- 100
  pars <- c(0.0001, 2.2, 0.005, 0.001, 1)
  area_pars <- create_area_pars(
    max_area = 5000,
    current_area = 2500,
    proportional_peak_t = 0.5,
    total_island_age = 15,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0)
  nonoceanic_pars <- c(0, 0)
  hyper_pars <- create_hyper_pars(d = 0.2, x = 0.1)
  island_ontogeny <- DAISIE::translate_island_ontogeny("beta")
  sea_level <- DAISIE:::translate_sea_level("const")
  peak <- DAISIE:::calc_peak(totaltime = totaltime,
                             area_pars = area_pars)
  Amax <- DAISIE:::get_global_max_area(totaltime = totaltime,
                                       area_pars = area_pars,
                                       peak = peak,
                                       island_ontogeny = island_ontogeny,
                                       sea_level = sea_level)
  Amin <- DAISIE:::get_global_min_area(totaltime = totaltime,
                                       area_pars = area_pars,
                                       peak = peak,
                                       island_ontogeny = island_ontogeny,
                                       sea_level = sea_level)
  set.seed(234567890)
  expect_silent(
    DAISIE:::DAISIE_sim_core_time_dependent(
      time = totaltime,
      mainland_n = mainland_n,
      pars = pars,
      area_pars = area_pars,
      hyper_pars = hyper_pars,
      nonoceanic_pars = nonoceanic_pars,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      peak = peak,
      Amax = Amax,
      Amin = Amin,
      extcutoff = 1000
    )
  )
})

test_that("Ontogeny oceanic should run silent CS", {
  totaltime <- 10
  mainland_n <- 1
  pars <- c(0.0001, 2.2, 0.005, 0.001, 1)
  island_ontogeny <- DAISIE::translate_island_ontogeny("beta")
  sea_level <- DAISIE:::translate_sea_level("const")
  area_pars <- create_area_pars(
    max_area = 5000,
    current_area = 2500,
    proportional_peak_t = 0.5,
    total_island_age = 15,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0)
  nonoceanic_pars <- c(0, 0)
  hyper_pars <- create_hyper_pars(d = 0.2, x = 0.1)
  peak <- DAISIE:::calc_peak(totaltime = totaltime,
                             area_pars = area_pars)
  Amax <- DAISIE:::get_global_max_area(totaltime = totaltime,
                                       area_pars = area_pars,
                                       peak = peak,
                                       island_ontogeny = island_ontogeny,
                                       sea_level = sea_level)
  Amin <- DAISIE:::get_global_min_area(totaltime = totaltime,
                                       area_pars = area_pars,
                                       peak = peak,
                                       island_ontogeny = island_ontogeny,
                                       sea_level = sea_level)
  set.seed(420)
  expect_silent(
    DAISIE:::DAISIE_sim_core_time_dependent(
      time = totaltime,
      mainland_n = mainland_n,
      pars = pars,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      area_pars = area_pars,
      hyper_pars = hyper_pars,
      nonoceanic_pars = nonoceanic_pars,
      peak = peak,
      Amax = Amax,
      Amin = Amin
    )
  )
})

test_that("Ontogeny oceanic with sea level should run silent CS", {
  totaltime <- 10
  mainland_n <- 1
  pars <- c(0.00001, 2.2, 0.005, 0.06, 1)
  area_pars <- create_area_pars(
    max_area = 5000,
    current_area = 2500,
    proportional_peak_t = 0.5,
    total_island_age = 15,
    sea_level_amplitude = 60,
    sea_level_frequency = 10,
    island_gradient_angle = 30)
  nonoceanic_pars <- c(0, 0)
  hyper_pars <- create_hyper_pars(d = 0.2, x = 0.1)
  island_ontogeny <- translate_island_ontogeny("beta")
  sea_level <- translate_sea_level("const")
  peak <- DAISIE:::calc_peak(totaltime = totaltime,
                             area_pars = area_pars)
  Amax <- DAISIE:::get_global_max_area(totaltime = totaltime,
                                       area_pars = area_pars,
                                       peak = peak,
                                       island_ontogeny = island_ontogeny,
                                       sea_level = sea_level)
  Amin <- DAISIE:::get_global_min_area(totaltime = totaltime,
                                       area_pars = area_pars,
                                       peak = peak,
                                       island_ontogeny = island_ontogeny,
                                       sea_level = sea_level)
  set.seed(439)
  expect_silent(
    out <- DAISIE:::DAISIE_sim_core_time_dependent(
      time = totaltime,
      mainland_n = mainland_n,
      pars = pars,
      area_pars = area_pars,
      hyper_pars = hyper_pars,
      nonoceanic_pars = nonoceanic_pars,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      peak = peak,
      Amax = Amax,
      Amin = Amin
    )
  )
  expect_equal(out$branching_times, c(10, 0.17840243993784999))
})

test_that("all species extinct if island dead", {
  totaltime <- 10
  mainland_n <- 1000
  pars <- c(0.0001, 100, 0.005, 0.0001, 0.1)
  hyper_pars <- create_hyper_pars(d = 0.2, x = 0.1)
  area_pars <- create_area_pars(
    max_area = 5000,
    current_area = 1,
    proportional_peak_t = 0.5,
    total_island_age = 10.1,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0)
  nonoceanic_pars <- c(0, 0)
  island_ontogeny <- DAISIE::translate_island_ontogeny("beta")
  sea_level <- DAISIE:::translate_sea_level("const")
  peak <- DAISIE:::calc_peak(totaltime = totaltime,
                             area_pars = area_pars)
  Amax <- DAISIE:::get_global_max_area(totaltime = totaltime,
                                       area_pars = area_pars,
                                       peak = peak,
                                       island_ontogeny = island_ontogeny,
                                       sea_level = sea_level)
  Amin <- DAISIE:::get_global_min_area(totaltime = totaltime,
                                       area_pars = area_pars,
                                       peak = peak,
                                       island_ontogeny = island_ontogeny,
                                       sea_level = sea_level)
  set.seed(1)
  ontogeny_sim <- DAISIE:::DAISIE_sim_core_time_dependent(
    time = totaltime,
    mainland_n = mainland_n,
    pars = pars,
    area_pars = area_pars,
    hyper_pars = hyper_pars,
    nonoceanic_pars = nonoceanic_pars,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    peak = peak,
    Amax = Amax,
    Amin = Amin,
    extcutoff = 1000
  )
  last_entry <- ontogeny_sim$stt_table[nrow(ontogeny_sim$stt_table), ]
  expect_true(last_entry[1] == 0)
  expect_true(last_entry[2] == 0)
  expect_true(last_entry[3] == 0)
  expect_true(last_entry[4] == 0)
})

test_that("!is.null(area_pars) && island_ontogeny == 'const'", {
  totaltime <- 1
  mainland_n <- 100
  pars <- c(2, 2, 20, 0, 1)
  island_ontogeny <- 0
  sea_level <- 0
  area_pars <- create_area_pars(
    max_area = 1,
    current_area = 1,
    proportional_peak_t = 1,
    total_island_age = 1,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0
  )
  nonoceanic_pars <- c(0, 0)
  hyper_pars <- create_hyper_pars(d = 0, x = 0)
  peak <- 1
  Amax <- 1
  Amin <- 1
  expect_error(
    DAISIE:::DAISIE_sim_core_time_dependent(
      time = totaltime,
      mainland_n = mainland_n,
      pars = pars,
      island_ontogeny = 0,
      sea_level = 0,
      peak = peak,
      Amin = Amin,
      Amax = Amax,
      area_pars = area_pars,
      hyper_pars = hyper_pars,
      nonoceanic_pars = nonoceanic_pars
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
                island_ontogeny = island_ontogeny,
                sea_level = sea_level
              ), regexp =
                "Island ontogeny and/or sea level specified but area parameters not
    available. Please either set island_ontogeny and sea_level to NULL, or
    specify area_pars"
            )
          })
