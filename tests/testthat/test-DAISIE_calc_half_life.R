context("DAISIE_calc_half_life")
test_that("DAISIE_calc_half_life gives a valid half life", {
  time <- 1
  mainland_n <- 10
  pars <- c(1, 1, 20, 0.1, 1)
  island_replicates <- DAISIE_sim_core(time = time,
                                       mainland_n = mainland_n,
                                       pars = pars)
  half_life <- DAISIE_calc_half_life(island_replicates, pars)
  expect_true(is.numeric(half_life))
  expect_gt(half_life, 0)
})

test_that("DAISIE_calc_half_life gives a vector of half lives
          for oceanic IW simulation replicates", {
  time <- 1
  mainland_n <- 10
  pars <- c(1, 1, 20, 0.1, 1)
  replicates <- 10
  island_replicates <- list()
  for (rep in 1:replicates) {
    island_replicates[[rep]] <- DAISIE_sim_core(time = time,
                                                mainland_n = mainland_n,
                                                pars = pars)
  }
    half_life <- DAISIE_calc_half_life(island_replicates, pars)
    expect_true(is.numeric(half_life))
                expect_true(is.vector(half_life))
                expect_length(half_life, replicates)
})

test_that("DAISIE_calc_half_life gives a vector of half lives
          for oceanic CS simulation replicates", {
  time <- 1
  mainland_n <- 10
  pars <- c(1, 1, 20, 0.1, 1)
  replicates <- 10
  for (rep in 1:replicates) {
    island_replicates[[rep]] <- list()
    full_list <- list()
    for (m_spec in 1:mainland_n) {
      full_list[[m_spec]] <- DAISIE_sim_core(
        time = time,
        mainland_n = 1,
        pars = pars)
      }
    island_replicates[[rep]] <- full_list
    }
  half_life <- DAISIE_calc_half_life(island_replicates, pars)
  expect_true(is.numeric(half_life))
  expect_true(is.vector(half_life))
  expect_length(half_life, replicates)
})

test_that("DAISIE_calc_half_life gives a vector of half lives
          for non-oceanic IW simulation replicates", {
  time <- 1
  mainland_n <- 10
  pars <- c(1, 1, 20, 0.1, 1)
  replicates <- 10
  island_replicates <- list()
  for (rep in 1:replicates) {
    island_replicates[[rep]] <- DAISIE_sim_core(time = time,
                                                mainland_n = mainland_n,
                                                pars = pars,
                                                island_type = "nonoceanic",
                                                nonoceanic = c(0.1, 0.9))
  }
  half_life <- DAISIE_calc_half_life(island_replicates, pars)
  expect_true(is.numeric(half_life))
  expect_true(is.vector(half_life))
  expect_length(half_life, replicates)
})

test_that("DAISIE_calc_half_life gives a vector of half lives
          for non-oceanic CS simulation replicates", {
  time <- 1
  mainland <- 10
  pars <- c(1, 1, 20, 0.1, 1)
  replicates <- 10
  for (rep in 1:replicates) {
    island_replicates[[rep]] <- list()
    full_list <- list()
    for (m_spec in 1:mainland_n) {
      full_list[[m_spec]] <- DAISIE_sim_core(
        time = time,
        mainland_n = 1,
        pars = pars,
        island_type = "nonoceanic",
        nonoceanic = c(0.1, 0.9))
    }
    island_replicates[[rep]] <- full_list
  }
  half_life <- DAISIE_calc_half_life(island_replicates, pars)
  expect_true(is.numeric(half_life))
  expect_true(is.vector(half_life))
  expect_length(half_life, replicates)
})

test_that("Half-life becomes shorter for diversity-dependent extinction", {
  time <- 1
  mainland_n <- 10
  pars <- c(1, 1, 20, 0.1, 1)
  ddmodel <- c(1, 1, 1)
  DD_sim <- DAISIE_sim_core(time = time,
                            mainland_n = mainland_n,
                            pars = pars,
                            ddmodel = ddmodel)
  time <- 1
  mainland_n <- 10
  pars <- c(1, 1, 20, 0.1, 1)
  ddmodel <- c(1, 0, 1)
  DI_sim <- DAISIE_sim_core(time = time,
                            mainland_n = mainland_n,
                            pars = pars,
                            ddmodel = ddmodel)
  DD_half_life <- DAISIE_calc_half_life(DD_sim, pars)
  DI_half_life <- DAISIE_calc_half_life(DI_sim, pars)
  expect_gt(DI_half_life, DD_half_life)
  expect_true(DI_half_life != DD_half_life)
})