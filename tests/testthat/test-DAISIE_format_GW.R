context("test-DAISIE_format_GW")

test_that("use GW with empty island", {
  skip("NEEDS FIXING ON BRANCH")
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
    num_guilds = num_guilds
  )
  island_replicates[[1]] <- out
  expect_silent(
    formated_GW_sim <- DAISIE:::DAISIE_format_GW(
      island_replicates = island_replicates,
      time = time,
      M = mainland_n,
      sample_freq = sample_freq,
      island_type = island_type,
      verbose = verbose,
      start_midway = start_midway
    )
  )
})

test_that("use mkGW with empty island", {
  skip("NEEDS FIXING ON BRANCH")
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
    num_guilds = num_guilds
  )
  island_replicates[[1]] <- out
  expect_silent(
    formated_GW_sim <- DAISIE:::DAISIE_format_GW(
      island_replicates = island_replicates,
      time = time,
      M = mainland_n,
      sample_freq = sample_freq,
      island_type = island_type,
      verbose = verbose,
      start_midway = start_midway
    )
  )
})


test_that("use GW with non-empty island", {
  skip("NEEDS FIXING ON BRANCH")
  pars <- c(0.4, 0.1, 10, 1, 0.5)
  time <- 1
  mainland_n <- 1
  island_type <- "oceanic"
  verbose <- FALSE
  sample_freq <- 1
  start_midway <- FALSE
  set.seed(1)
  island_replicates <- list()
  out <- list()
  out[[1]] <- DAISIE:::DAISIE_sim_core(
    time = time,
    pars = pars,
    mainland_n
    = mainland_n
  )
  island_replicates[[1]] <- out
  expect_silent(
    formated_CS_sim <- DAISIE:::DAISIE_format_CS(
      island_replicates = island_replicates,
      time = time,
      M = mainland_n,
      sample_freq = sample_freq,
      island_type = island_type,
      verbose = verbose,
      start_midway = start_midway
    )
  )
})


test_that("use mkGW with non-empty island", {
  skip("NEEDS FIXING ON BRANCH")
  pars <- c(0.4, 0.1, 10, 1, 0.5)
  time <- 1
  mainland_n <- 1
  island_type <- "oceanic"
  verbose <- FALSE
  sample_freq <- 1
  start_midway <- FALSE
  set.seed(1)
  island_replicates <- list()
  out <- list()
  out[[1]] <- DAISIE:::DAISIE_sim_core(
    time = time,
    pars = pars,
    mainland_n
    = mainland_n
  )
  island_replicates[[1]] <- out
  expect_silent(
    formated_CS_sim <- DAISIE:::DAISIE_format_CS(
      island_replicates = island_replicates,
      time = time,
      M = mainland_n,
      sample_freq = sample_freq,
      island_type = island_type,
      verbose = verbose,
      start_midway = start_midway
    )
  )
})

test_that("abuse", {
  skip("NEEDS FIXING ON BRANCH")
  expect_error(DAISIE:::DAISIE_format_GW("nonsense"))
})

