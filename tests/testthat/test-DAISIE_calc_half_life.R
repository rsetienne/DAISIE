context("DAISIE_calc_half_life")

test_that("DAISIE_calc_half_life gives a vector of half lives
          for oceanic IW simulation replicates", {
  time <- 2
  mainland_n <- 100
  pars <- c(2, 1, 40, 0.1, 1)
  replicates <- 1
  divdepmodel <- "IW"
  island_replicates <- list()
  set.seed(1)
  for (rep in 1:replicates) {
    island_replicates[[rep]] <- DAISIE_sim_core(time = time,
                                                mainland_n = mainland_n,
                                                pars = pars)
    
  }
    half_life <- DAISIE_calc_half_life(island_replicates,
                                       mainland_n,
                                       pars,
                                       island_type,
                                       divdepmodel)
    expect_true(is.list(half_life))
    expect_true(is.numeric(half_life[[1]]))
    expect_gt(half_life[[1]], 0)
    expect_length(half_life, replicates)
})

test_that("DAISIE_calc_half_life gives a vector of half lives
          for oceanic CS simulation replicates", {

  # Test takes too long, only run on Travis              
  if (Sys.getenv("TRAVIS") != "") return()

  time <- 2
  mainland_n <- 10
  pars <- c(3, 1, 40, 0.5, 1)
  replicates <- 1
  divdepmodel <- "CS"
  island_replicates <- list()
  set.seed(1)
  for (rep in 1:replicates) {
    island_replicates[[rep]] <- list()
    full_list <- list()
    for (m_spec in 1:mainland_n) {
      full_list[[m_spec]] <- DAISIE_sim_core(time = time,
                                             mainland_n = 1,
                                             pars = pars)
      }
    island_replicates[[rep]] <- full_list
    }
  half_life <- DAISIE_calc_half_life(island_replicates,
                                     mainland_n,
                                     pars,
                                     divdepmodel = divdepmodel)
  expect_true(is.list(half_life))
  expect_true(is.numeric(half_life[[1]]))
  expect_gt(half_life[[1]], 0)
  expect_length(half_life, replicates)
})

test_that("DAISIE_calc_half_life gives a vector of half lives
          for non-oceanic IW simulation replicates", {
  time <- 2
  mainland_n <- 1000
  pars <- c(2, 2, 20, 0.1, 1)
  replicates <- 1
  divdepmodel <- "IW"
  island_replicates <- list()
  set.seed(1)
  for (rep in 1:replicates) {
    island_replicates[[rep]] <- DAISIE_sim_core(time = time,
                                                mainland_n = mainland_n,
                                                pars = pars,
                                                island_type = "nonoceanic",
                                                nonoceanic = c(0.1, 0.9))
  }
  half_life <- DAISIE_calc_half_life(island_replicates,
                                     mainland_n,
                                     pars,
                                     island_type = "nonoceanic",
                                     divdepmodel)
  expect_true(is.list(half_life))
  expect_true(is.numeric(half_life[[1]]))
  expect_gt(half_life[[1]], 0)
  expect_length(half_life, replicates)
})

test_that("DAISIE_calc_half_life gives a vector of half lives
          for non-oceanic CS simulation replicates", {
skip("test takes too long")
  time <- 2
  mainland <- 10
  pars <- c(2, 2, 2, 0.1, 1)
  replicates <- 1
  divdepmodel <- "CS"
  island_replicates <- list()
  set.seed(1)
  for (rep in 1:replicates) {
    island_replicates[[rep]] <- list()
    full_list <- list()
    for (m_spec in 1:mainland_n) {
      full_list[[m_spec]] <- DAISIE_sim_core(
        time = time,
        mainland_n = 1,
        pars = pars,
        island_type = "nonoceanic",
        nonoceanic = c(0.9, 0.9))
    }
    island_replicates[[rep]] <- full_list
  }
  half_life <- DAISIE_calc_half_life(island_replicates,
                                     mainland_n,
                                     pars,
                                     island_type = "nonoceanic",
                                     divdepmodel)
  expect_true(is.list(half_life))
  expect_true(is.numeric(half_life[[1]]))
  expect_gt(half_life[[1]], 0)
  expect_length(half_life, replicates)
})

test_that("Half-life becomes shorter for diversity-dependent extinction", {
skip("needs fixing on branch")
  time <- 1
  mainland_n <- 10
  pars <- c(1, 1, 20, 0.1, 1)
  ddmodel_sim <- 11
  divdepmodel <- "IW"
  DD_sim <- DAISIE_sim_core(time = time,
                            mainland_n = mainland_n,
                            pars = pars,
                            ddmodel_sim = ddmodel_sim)
  time <- 1
  mainland_n <- 10
  pars <- c(1, 1, 20, 0.1, 1)
  ddmodel_sim <- 11
  divdepmodel <- "IW"
  DI_sim <- DAISIE_sim_core(time = time,
                            mainland_n = mainland_n,
                            pars = pars,
                            ddmodel_sim = ddmodel_sim)
  DD_half_life <- DAISIE_calc_half_life(DD_sim, pars, divdepmodel)
  DI_half_life <- DAISIE_calc_half_life(DI_sim, pars, divdepmodel)
  expect_gt(DI_half_life, DD_half_life)
  expect_true(DI_half_life != DD_half_life)
})

test_that("Expect NA when extinction is zero for nonoceanic islands as half-life
          will not be reached", {
  time <- 1
  mainland_n <- 10
  pars <- c(1, 0, 20, 0.1, 1)
  replicates <- 1
  divdepmodel <- "IW"
  island_replicates <- list()
  set.seed(1)
  for (rep in 1:replicates) {
    island_replicates[[rep]] <- DAISIE_sim_core(time = time,
                                                mainland_n = mainland_n,
                                                pars = pars,
                                                island_type = "nonoceanic",
                                                nonoceanic = c(0.1, 0.9))
  }
  half_life <- DAISIE_calc_half_life(island_replicates,
                                     mainland_n,
                                     pars,
                                     island_type,
                                     divdepmodel)
  expect_identical(half_life[[1]][[1]], NA)
})
