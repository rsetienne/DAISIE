context("test-DAISIE_format_IW")

test_that("use with empty island", {
  skip("NEEDS FIXING ON BRANCH")
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
})

test_that("use with non-empty island", {
  skip("NEEDS FIXING ON BRANCH")
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
  skip("NEEDS FIXING ON BRANCH")
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
  skip("NEEDS FIXING ON BRANCH")
  expect_error(DAISIE:::DAISIE_format_IW(
    island_replicates = "nonsense",
    time = time,
    M = mainland_n,
    sample_freq = sample_freq,
    verbose = verbose
  ))
})
