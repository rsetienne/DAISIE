context("DAISIE_sim_constant_rate_shift")

test_that("use split-rates model", {
  expect_silent(
    DAISIE_sim_constant_rate_shift(
      time = 10,
      M = 10,
      pars = c(1, 1, 1, 0.1, 1, 1, 1, 1, 0.1, 1),
      replicates = 1,
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

