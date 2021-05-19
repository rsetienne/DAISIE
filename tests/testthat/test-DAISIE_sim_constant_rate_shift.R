context("DAISIE_sim_constant_rate_shift")

test_that("use CS split-rates model", {
  expect_silent(DAISIE_sim_constant_rate_shift(
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

test_that("us CS split-rates with cond", {
  set.seed(Sys.time()) # Always run a different sim
  time <- 10
  M <- 10
  pars <- c(1, 1, 1, 0.1, 1, 1, 1, 1, 0.1, 1)
  replicates <- 1
  divdepmodel <- "CS"
  shift_times <- 5
  cond <- 5
  expect_silent(
    out <- DAISIE_sim_constant_rate_shift(
      time = time,
      M = M,
      pars = pars,
      replicates = 1,
      divdepmodel = "CS",
      shift_times = shift_times,
      plot_sims = FALSE,
      verbose = FALSE,
      cond = cond
    )
  )

  expect_true(out[[1]][[1]]$stt_all[nrow(out[[1]][[1]]$stt_all), 5] >= cond)
})

test_that("expected cond or 0 cond CS split-rates model", {
  set.seed(1)
  cond <- 5
  expect_silent(out_no_cond <- DAISIE_sim_constant_rate_shift(
    time = 10,
    M = 10,
    pars = c(1, 1, 1, 0.1, 1, 1, 1, 1, 0.1, 1),
    replicates = 1,
    divdepmodel = "CS",
    shift_times = 5,
    plot_sims = FALSE,
    verbose = FALSE
  ))

  expect_true(
    out_no_cond[[1]][[1]]$stt_all[nrow(out_no_cond[[1]][[1]]$stt_all), 5] < cond
  )

  set.seed(1)
  expect_silent(out_cond <- DAISIE_sim_constant_rate_shift(
    time = 10,
    M = 10,
    pars = c(1, 1, 1, 0.1, 1, 1, 1, 1, 0.1, 1),
    replicates = 1,
    divdepmodel = "CS",
    shift_times = 5,
    plot_sims = FALSE,
    verbose = FALSE,
    cond = cond
  ))

  expect_true(
    out_cond[[1]][[1]]$stt_all[nrow(out_cond[[1]][[1]]$stt_all), 5] >= cond
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


test_that("Reference output matches DAISIE_sim_constant_rate_shift ", {
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
    unname(sims[[1]][[1]]$stt_all[26, ]), c(0, 42, 9, 2, 51)
  )
})

test_that("The SR simulation code works", {

  # Simulate fish diversity over 4 Ma
  set.seed(1)
  M <- 312
  IslandAge <- 4
  sims <- DAISIE_SR_sim(
    time = 4,
    M = M - 17,
    pars = pars1,
    replicates = 1,
    plot_sims = FALSE,
    ddep = 11
  )
  # Compare richnesses of the last time bin
  testthat::expect_equal(
    unname(sims[[1]][[1]]$stt_all[26, ]), c(0, 56, 11, 0, 66)
  )
})
