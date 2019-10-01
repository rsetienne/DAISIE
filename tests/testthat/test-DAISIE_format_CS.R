context("test-DAISIE_format_CS")

test_that("use with empty island", {
  pars <- c(0.4, 0.2, 10, 0.0001, 0.5)
  time <- 1
  mainland_n <- 1
  verbose <- FALSE
  sample_freq <- 1
  start_midway <- FALSE
  Tpars <- NULL
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
    formated_CS_sim <- DAISIE:::DAISIE_format_CS(
    island_replicates = island_replicates,
    time = time,
    M = mainland_n,
    Tpars = Tpars,
    sample_freq = sample_freq,
    verbose = verbose,
    start_midway = start_midway
    )
  )
})

test_that("use with non-empty island", {
  pars <- c(0.4, 0.1, 10, 1, 0.5)
  time <- 1
  mainland_n <- 1
  verbose <- FALSE
  sample_freq <- 1
  start_midway <- FALSE
  Tpars <- NULL
  set.seed(1)
  island_replicates <- list()
  out <- list()
  out[[1]] <- DAISIE:::DAISIE_sim_core(
    time = time,
    pars = pars,
    mainland_n = mainland_n,
    Tpars = Tpars
  )
  island_replicates[[1]] <- out

  expect_silent(
    formated_CS_sim <- DAISIE:::DAISIE_format_CS(
      island_replicates = island_replicates,
      time = time,
      M = mainland_n,
      Tpars = Tpars,
      sample_freq = sample_freq,
      verbose = verbose,
      start_midway = start_midway
    )
  )
})

test_that("abuse", {
  expect_error(DAISIE:::DAISIE_format_CS("nonsense"))
})
