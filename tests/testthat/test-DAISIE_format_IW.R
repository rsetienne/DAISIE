context("test-DAISIE_format_IW")

test_that("silent with empty island with correct output", {
  pars <- c(0.4, 0.2, 10, 0.0001, 0.5)
  time <- 1
  mainland_n <- 1000
  verbose <- FALSE
  sample_freq <- 1
  start_midway <- FALSE
  set.seed(1)
  island_replicates <- list()
  island_replicates[[1]] <- DAISIE:::DAISIE_sim_core(
    time = time,
    pars = pars,
    mainland_n = mainland_n
  )
  expect_silent(
    formated_IW_sim <- DAISIE:::DAISIE_format_IW(
      island_replicates = island_replicates,
      time = time,
      M = mainland_n,
      sample_freq = sample_freq,
      verbose = verbose,
      island_type = "oceanic"
    )
  )
  expected_IW_format <- list()
  expected_IW_format[[1]] <- list()
  stt_all <- matrix(ncol = 4, nrow = 2)
  colnames(stt_all) <- c("Time", "nI", "nA", "nC")
  stt_all[1, ] <- c(1, 0, 0, 0)
  stt_all[2, ] <- c(0, 0, 0, 0)
  brts_table <- matrix(ncol = 4, nrow = 1)
  colnames(brts_table) <- c("brt", "clade", "event", "endemic")
  brts_table[1, ] <- c(10, 0, 0, NA)
  expected_IW_format[[1]][[1]] <- list(island_age = 1,
                                       not_present = 1000,
                                       stt_all = stt_all,
                                       init_nonend_spec = 0,
                                       init_end_spec = 0,
                                       brts_table = brts_table)
  expect_equal(formated_IW_sim, expected_IW_format)
})

test_that("silent with non-empty island with correct output", {
  pars <- c(0.4, 0.2, 10, 1, 0.5)
  time <- 1
  mainland_n <- 1000
  verbose <- FALSE
  sample_freq <- 1
  start_midway <- FALSE
  set.seed(1)
  island_replicates <- list()
  island_replicates[[1]] <- DAISIE:::DAISIE_sim_core(
    time = time,
    pars = pars,
    mainland_n = mainland_n
  )
  expect_silent(
    formated_IW_sim <- DAISIE:::DAISIE_format_IW(
      island_replicates = island_replicates,
      time = time,
      M = mainland_n,
      sample_freq = sample_freq,
      verbose = verbose,
      island_type = "oceanic"
    )
  )
})

test_that("DAISIE_format_IW prints when verbose = TRUE", {
  pars <- c(0.4, 0.2, 10, 1, 0.5)
  time <- 1
  mainland_n <- 1000
  verbose <- TRUE
  sample_freq <- 1
  start_midway <- FALSE
  set.seed(1)
  island_replicates <- list()
  island_replicates[[1]] <- DAISIE:::DAISIE_sim_core(
    time = time,
    pars = pars,
    mainland_n = mainland_n
  )
  expect_output(
    formated_IW_sim <- DAISIE:::DAISIE_format_IW(
      island_replicates = island_replicates,
      time = time,
      M = mainland_n,
      sample_freq = sample_freq,
      verbose = verbose,
      island_type = "oceanic"
    ),
    "Island being formatted: 1/1"
  )
})

test_that("Add_brt_table [insert verb] if (length(stac1_5s) != 0)", {
  skip("WIP")
  pars <- c(0.4, 0.2, 10, 1, 0.5)
  time <- 1
  mainland_n <- 1000
  verbose <- FALSE
  sample_freq <- 1
  start_midway <- FALSE
  set.seed(1)
  island_replicates <- list()
  island_replicates[[1]] <- DAISIE:::DAISIE_sim_core(
    time = time,
    pars = pars,
    mainland_n = mainland_n
  )
  expect_(
    formated_IW_sim <- DAISIE:::DAISIE_format_IW(
      island_replicates = island_replicates,
      time = time,
      M = mainland_n,
      sample_freq = sample_freq,
      verbose = verbose,
      island_type = "oceanic",
      verbose = verbose
    )
  )
})

test_that("abuse", {
  expect_error(DAISIE:::DAISIE_format_IW("nonsense"))
})
