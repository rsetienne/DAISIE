context("DAISIE_sim_constant_rate_shift")

test_that("use CS split-rates model", {
  expect_silent(
    DAISIE_sim_constant_rate_shift(
      time = 10,
      M = 10,
      pars = c(1, 1, 1, 0.1, 1, 1, 1, 1, 0.1, 1),
      replicates = 1,
      divdepmodel = "CS",
      shift_times = 5,
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
})

test_that("use IW split-rates model", {
  expect_silent(
    DAISIE_sim_constant_rate_shift(
      time = 10,
      M = 10,
      pars = c(1, 1, 1, 0.1, 1, 1, 1, 1, 0.1, 1),
      replicates = 1,
      divdepmodel = "IW",
      shift_times = 5,
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
})

test_that("use GW split-rates model", {
  expect_silent(
    DAISIE_sim_constant_rate_shift(
      time = 10,
      M = 10,
      pars = c(1, 1, 1, 0.1, 1, 1, 1, 1, 0.1, 1),
      replicates = 1,
      divdepmodel = "GW",
      num_guilds = 2,
      shift_times = 5,
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
})

test_that("abuse split-rates model", {
  expect_error(DAISIE_sim_constant_rate_shift(
    time = 1,
    M = 1,
    pars = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    replicates = 1,
    shift_times = 5,
    verbose = FALSE,
    plot_sims = FALSE
  ))
  expect_error(DAISIE_sim_constant_rate_shift(
    time = 10,
    M = 1,
    pars = c(1, 1, 1, 1, 1),
    replicates = 1,
    shift_times = 5,
    verbose = FALSE,
    plot_sims = FALSE
  ))
})

test_that("split-rates model prints when verbose = TRUE", {
  expect_output(
    DAISIE_sim_constant_rate_shift(
      time = 10,
      M = 10,
      pars = c(1, 1, 1, 0.1, 1, 1, 1, 1, 0.1, 1),
      replicates = 1,
      shift_times = 5,
      plot_sims = FALSE,
      verbose = TRUE
    ),
    regexp = "Island replicate 1"
  )
})


test_that("testing the split_rate model is the same as before", {
  set.seed(1)
  M <- 312
  island_age <- 4
  pars1 <- c(0.077, 0.956, Inf, 0.138, 0.442,
             0.077, 0.956, Inf, 0.655, 0.442)
  sims <- DAISIE_sim_constant_rate_shift(
    time = island_age,
    M = 295,
    pars = pars1,
    replicates = 1,
    plot_sims = FALSE,
    shift_times = 0.1951,
    verbose = FALSE
  )
  # Compare richnesses of the last time bin
  testthat::expect_equal(
    unname(sims[[1]][[1]]$stt_all[26, ]), c(0, 56, 11, 0, 66)
  )
})

