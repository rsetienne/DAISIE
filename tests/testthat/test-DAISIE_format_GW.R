context("DAISIE_format_GW")

test_that("silent with empty island with correct output", {
  pars <- c(0.4, 0.2, 10, 0.0001, 0.5)
  time <- 1
  mainland_n <- 1
  num_guilds <- 1
  verbose <- FALSE
  sample_freq <- 1
  area_pars <- DAISIE::create_area_pars(
    max_area = 1,
    current_area = 1,
    proportional_peak_t = 0,
    total_island_age = 0,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0)
  hyper_pars <- create_hyper_pars(d = 0, x = 0)
  nonoceanic_pars <- c(0, 0)
  set.seed(1)
  island_replicates <- list()
  out <- list()
  out[[1]] <- DAISIE:::DAISIE_sim_core_cr(
    time = time,
    pars = pars,
    mainland_n = mainland_n,
    area_pars = area_pars,
    hyper_pars = hyper_pars,
    nonoceanic_pars = nonoceanic_pars
  )
  island_replicates[[1]] <- out
  expect_silent(
    formated_GW_sim <- DAISIE:::DAISIE_format_GW(
      island_replicates = island_replicates,
      time = time,
      M = mainland_n,
      sample_freq = sample_freq,
      num_guilds = num_guilds,
      verbose = verbose
    )
  )
  expected_GW_format <- list()
  expected_GW_format[[1]] <- list()
  stt_all <- matrix(ncol = 5, nrow = 2)
  colnames(stt_all) <- c("Time", "nI", "nA", "nC", "present")
  stt_all[1, ] <- c(1, 0, 0, 0, 0)
  stt_all[2, ] <- c(0, 0, 0, 0, 0)
  expected_GW_format[[1]][[1]] <- list(island_age = 1,
                                       not_present = 1,
                                       stt_all = stt_all)
  expect_identical(formated_GW_sim, expected_GW_format)
})

test_that("silent with non-empty island with correct output", {
  pars <- c(0.4, 0.1, 10, 1, 0.5)
  time <- 1
  mainland_n <- 1
  num_guilds <- 1
  verbose <- FALSE
  sample_freq <- 1
  area_pars <- DAISIE::create_area_pars(
    max_area = 1,
    current_area = 1,
    proportional_peak_t = 0,
    total_island_age = 0,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0)
  hyper_pars <- create_hyper_pars(d = 0, x = 0)
  nonoceanic_pars <- c(0, 0)
  set.seed(1)
  island_replicates <- list()
  out <- list()
  out[[1]] <- DAISIE:::DAISIE_sim_core_cr(
    time = time,
    pars = pars,
    mainland_n = mainland_n,
    area_pars = area_pars,
    hyper_pars = hyper_pars,
    nonoceanic_pars = nonoceanic_pars
  )
  island_replicates[[1]] <- out
  expect_silent(
    formated_GW_sim <- DAISIE:::DAISIE_format_GW(
      island_replicates = island_replicates,
      time = time,
      M = mainland_n,
      num_guilds = num_guilds,
      sample_freq = sample_freq,
      verbose = verbose
    )
  )
})

test_that("output with empty island and verbose = TRUE", {
  pars <- c(0.4, 0.2, 10, 0.0001, 0.5)
  time <- 1
  mainland_n <- 1
  num_guilds <- 1
  verbose <- TRUE
  sample_freq <- 1
  set.seed(1)
  area_pars <- DAISIE::create_area_pars(
    max_area = 1,
    current_area = 1,
    proportional_peak_t = 0,
    total_island_age = 0,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0)
  hyper_pars <- create_hyper_pars(d = 0, x = 0)
  nonoceanic_pars <- c(0, 0)
  island_replicates <- list()
  out <- list()
  out[[1]] <- DAISIE:::DAISIE_sim_core_cr(
    time = time,
    pars = pars,
    mainland_n = mainland_n,
    area_pars = area_pars,
    hyper_pars = hyper_pars,
    nonoceanic_pars = nonoceanic_pars
  )
  island_replicates[[1]] <- out
  expect_output(
    formated_GW_sim <- DAISIE:::DAISIE_format_GW(
      island_replicates = island_replicates,
      time = time,
      M = mainland_n,
      sample_freq = sample_freq,
      num_guilds = num_guilds,
      verbose = verbose
    ),
    regexp = "Island being formatted: 1/1"
  )
})

test_that("abuse", {
  expect_error(DAISIE:::DAISIE_format_GW("nonsense"))
})
