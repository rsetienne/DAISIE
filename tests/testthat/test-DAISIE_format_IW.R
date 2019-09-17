context("test-DAISIE_format_IW")

test_that("use with empty island", {
  pars <- c(0.4, 0.2, 10, 0.0001, 0.5)
  time <- 1
  mainland_n <- 1000
  verbose <- FALSE
  sample_freq <- 1
  start_midway <- FALSE
  island_type <- "oceanic"
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
      Tpars = NULL,
      island_type = island_type
    )
  )
})

test_that("use with non-empty island", {
  pars <- c(0.4, 0.2, 10, 1, 0.5)
  time <- 1
  mainland_n <- 1000
  verbose <- FALSE
  sample_freq <- 1
  start_midway <- FALSE
  island_type <- "oceanic"
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
      Tpars = NULL,
      island_type = island_type
    )
  )
})

test_that("abuse", {
  expect_error(DAISIE:::DAISIE_format_IW(
    island_replicates = "nonsense",
    time = time,
    M = mainland_n,
    sample_freq = sample_freq,
    verbose = verbose,
    Tpars = NULL,
    island_type = "oceanic"
  ))
})
