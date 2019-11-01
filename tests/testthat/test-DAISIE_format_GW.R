context("test-DAISIE_format_GW")

test_that("silent with empty island with correct output", {
  pars <- c(0.4, 0.2, 10, 0.0001, 0.5)
  time <- 1
  mainland_n <- 1
  island_type <- "oceanic"
  num_guilds <- 1
  verbose <- FALSE
  sample_freq <- 1
  start_midway <- FALSE
  set.seed(1)
  island_replicates <- list()
  out <- list()
  out[[1]] <- DAISIE:::DAISIE_sim_core(
    time = time,
    pars = pars,
    mainland_n = mainland_n,
  )
  island_replicates[[1]] <- out
  expect_silent(
    formated_GW_sim <- DAISIE:::DAISIE_format_GW(
      island_replicates = island_replicates,
      time = time,
      M = mainland_n,
      sample_freq = sample_freq,
      island_type = island_type,
      num_guilds = num_guilds,
      start_midway = start_midway,
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
  expected_GW_format[[1]][[2]] <- list(branching_times = 1,
                                       stac = 0,
                                       missing_species = 0,
                                       init_nonend_spec = 0,
                                       init_end_spec = 0,
                                       carrying_capacity = "N/A",
                                       all_carrying_capacities = 10)
  expect_identical(formated_GW_sim, expected_GW_format)
})

test_that("multi-K silent with empty island with correct output", {
  pars <- c(0.4, 0.2, 10, 0.0001, 0.5)
  time <- 1
  mainland_n <- 1
  island_type <- "oceanic"
  num_guilds <- 1
  k_dist_params <- c(2, 0.5)
  verbose <- FALSE
  sample_freq <- 1
  start_midway <- FALSE
  set.seed(1)
  island_replicates <- list()
  out <- list()
  out[[1]] <- DAISIE:::DAISIE_sim_core(
    time = time,
    pars = pars,
    mainland_n = mainland_n,
    k_dist_params = k_dist_params
  )
  island_replicates[[1]] <- out
  expect_silent(
    formated_GW_sim <- DAISIE:::DAISIE_format_GW(
      island_replicates = island_replicates,
      time = time,
      M = mainland_n,
      num_guilds = num_guilds,
      sample_freq = sample_freq,
      island_type = island_type,
      verbose = verbose,
      start_midway = start_midway
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
  expected_GW_format[[1]][[2]] <- list(branching_times = 1,
                                       stac = 0,
                                       missing_species = 0,
                                       init_nonend_spec = 0,
                                       init_end_spec = 0,
                                       carrying_capacity = "N/A",
                                       all_carrying_capacities = 1.66173)
  expect_equal(formated_GW_sim, expected_GW_format)
})

test_that("silent with non-empty island with correct output", {
  pars <- c(0.4, 0.1, 10, 1, 0.5)
  time <- 1
  mainland_n <- 1
  island_type <- "oceanic"
  num_guilds <- 1
  verbose <- FALSE
  sample_freq <- 1
  start_midway <- FALSE
  set.seed(1)
  island_replicates <- list()
  out <- list()
  out[[1]] <- DAISIE:::DAISIE_sim_core(
    time = time,
    pars = pars,
    mainland_n = mainland_n
  )
  island_replicates[[1]] <- out
  expect_silent(
    formated_GW_sim <- DAISIE:::DAISIE_format_GW(
      island_replicates = island_replicates,
      time = time,
      M = mainland_n,
      num_guilds = num_guilds,
      sample_freq = sample_freq,
      island_type = island_type,
      verbose = verbose,
      start_midway = start_midway
    )
  )
})

test_that("multi-K silent with non-empty island with correct output", {
  pars <- c(0.4, 0.1, 10, 1, 0.5)
  time <- 1
  mainland_n <- 1
  island_type <- "oceanic"
  num_guilds <- 1
  k_dist_params <- c(2, 0.5)
  verbose <- FALSE
  sample_freq <- 1
  start_midway <- FALSE
  set.seed(1)
  island_replicates <- list()
  out <- list()
  out[[1]] <- DAISIE:::DAISIE_sim_core(
    time = time,
    pars = pars,
    mainland_n = mainland_n,
    k_dist_params = k_dist_params
  )
  island_replicates[[1]] <- out
  expect_silent(
    formated_GW_sim <- DAISIE:::DAISIE_format_GW(
      island_replicates = island_replicates,
      time = time,
      M = mainland_n,
      sample_freq = sample_freq,
      island_type = island_type,
      num_guilds = num_guilds,
      verbose = verbose,
      start_midway = start_midway
    )
  )
})


test_that("output with empty island and verbose = TRUE", {
  pars <- c(0.4, 0.2, 10, 0.0001, 0.5)
  time <- 1
  mainland_n <- 1
  island_type <- "oceanic"
  num_guilds <- 1
  verbose <- TRUE
  sample_freq <- 1
  start_midway <- FALSE
  set.seed(1)
  island_replicates <- list()
  out <- list()
  out[[1]] <- DAISIE:::DAISIE_sim_core(
    time = time,
    pars = pars,
    mainland_n = mainland_n,
  )
  island_replicates[[1]] <- out
  expect_output(
    formated_GW_sim <- DAISIE:::DAISIE_format_GW(
      island_replicates = island_replicates,
      time = time,
      M = mainland_n,
      sample_freq = sample_freq,
      island_type = island_type,
      num_guilds = num_guilds,
      start_midway = start_midway,
      verbose = verbose
    )
  )
})

test_that("abuse", {
  expect_error(DAISIE:::DAISIE_format_GW("nonsense"))
})
