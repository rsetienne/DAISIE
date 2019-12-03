context("DAISIE_sim split rate and DAISIE_SR_sim comparison")

test_that("testing the split_rate model is the same as before", {
  skip("WIP")
  set.seed(1)
  M <- 312
  island_age <- 4
  pars1 <- c(0.077, 0.956, Inf, 0.138, 0.442,
             0.077, 0.956, Inf, 0.655, 0.442)
  sims <- DAISIE_sim(
      time = island_age,
      M = 295,
      pars = pars1,
      replicates = 1,
      plot_sims = FALSE,
      shift_times = 0.1951,
      verbose = FALSE
    )


  #old code
  set.seed(1)
  M <- 312
  IslandAge <- 4
  sims <- DAISIE_sim(
    time = 4,
    M = M - 17,
    pars = pars1[1:10],
    replicates = 1,
    plot_sims = FALSE,
    pars_shift = TRUE,
    shift_times = 0.1951,
    verbose = FALSE
  )
  # Compare richnesses of the last time bin
  testthat::expect_equal(
    unname(sims[[1]][[1]]$stt_all[26, ]), c(0, 56, 11, 0, 66)
  )

})

