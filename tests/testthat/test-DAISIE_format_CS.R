context("test-DAISIE_format_CS")

test_that("silent with empty island with correct output", {
  pars <- c(0.4, 0.2, 10, 0.0001, 0.5)
  time <- 1
  mainland_n <- 1
  verbose <- FALSE
  sample_freq <- 1
  island_type = "oceanic"
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
    sample_freq = sample_freq,
    island_type = island_type,
    verbose = verbose
    )
  )
  expected_CS_format <- list()
  expected_CS_format[[1]] <- list()
  stt_all <- matrix(ncol = 5, nrow = 2)
  colnames(stt_all) <- c("Time", "nI", "nA", "nC", "present")
  stt_all[1, ] <- c(1, 0, 0, 0, 0)
  stt_all[2, ] <- c(0, 0, 0, 0, 0)
  expected_CS_format[[1]][[1]] <- list(island_age = 1,
                                  not_present = 1,
                                  stt_all = stt_all)
  expected_CS_format[[1]][[2]] <- list(branching_times = 1,
                                  stac = 0,
                                  missing_species = 0,
                                  init_nonend_spec = 0,
                                  init_end_spec = 0,
                                  carrying_capacity = "N/A",
                                  all_carrying_capacities = 10)
  expect_identical(formated_CS_sim, expected_CS_format)
})

test_that("silent with non-empty island with correct output", {
  pars <- c(0.5, 0.1, 10, 1, 0.5)
  time <- 1
  mainland_n <- 1
  verbose <- FALSE
  sample_freq <- 1
  island_type <- "oceanic"
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
      sample_freq = sample_freq,
      island_type = island_type,
      verbose = verbose
    )
  )
  expected_CS_format <- list()
  expected_CS_format[[1]] <- list()
  stt_all <- matrix(ncol = 5, nrow = 2)
  colnames(stt_all) <- c("Time", "nI", "nA", "nC", "present")
  stt_all[1, ] <- c(1, 0, 0, 0, 0)
  stt_all[2, ] <- c(0, 0, 0, 3, 1)
  expected_CS_format[[1]][[1]] <- list(island_age = 1,
                                       not_present = 0,
                                       stt_all = stt_all)
  expected_CS_format[[1]][[2]] <- list(branching_times = c(1.00000000,
                                                           0.24481817,
                                                           0.17312829,
                                                           0.02966824),
                                       stac = 2,
                                       missing_species = 0,
                                       init_nonend_spec = 0,
                                       init_end_spec = 0,
                                       carrying_capacity = 10)
  expected_CS_format[[1]][[3]] <- list(init_nonend_spec = 0,
                                       init_end_spec = 0,
                                       all_carrying_capacities = 10)
  expect_equal(formated_CS_sim, expected_CS_format)
})

test_that("output with empty island and verbose = TRUE", {
  pars <- c(0, 1, 1, 0.0001, 0)
  time <- 1
  mainland_n <- 1
  island_type <- "oceanic"
  verbose <- TRUE
  sample_freq <- 1
  set.seed(1)
  island_replicates <- list()
  out <- list()
  out[[1]] <- DAISIE:::DAISIE_sim_core(
    time = time,
    pars = pars,
    mainland_n = mainland_n
  )
  island_replicates[[1]] <- out
  expect_output(
    formated_CS_sim <- DAISIE:::DAISIE_format_CS(
      island_replicates = island_replicates,
      time = time,
      M = mainland_n,
      sample_freq = sample_freq,
      island_type = island_type,
      verbose = verbose
    )
  )
})

test_that("silent with empty 2 type island", {
  skip("NEEDS FINISHING ON BRANCH")
  pars <- c(0, 3, 1, 0.001, 0, 0, 3, 1, 0.001, 0)
  time <- 1
  M <- 1
  replicates <- 1
  prop_type2_pool <- 0.1
  replicates_apply_type2 <- TRUE
  verbose <- FALSE
  sample_freq <- 1
  set.seed(1)
  island_replicates <- list()
  out <- list()
  out[[1]] <- DAISIE:::DAISIE_sim_min_type2(
    time = totaltime,
    M = M,
    pars = pars,
    replicates = replicates,
    prop_type2_pool = prop_type2_pool,
    verbose = verbose)
  island_replicates[[1]] <- out
  expect_output(
    formated_CS_sim <- DAISIE:::DAISIE_format_CS(
      island_replicates = island_replicates,
      time = time,
      M = mainland_n,
      sample_freq = sample_freq,
      island_type = island_type,
      verbose = verbose
    )
  )
})

test_that("silent with non-empty 2 type island", {
  skip("NEEDS FINISHING ON BRANCH")
  pars <- c(0.4, 0.1, 10, 1, 0.5, 0.4, 0.1, 10, 1, 0.5)
  time <- 1
  M <- 10
  island_type <- "oceanic"
  verbose <- FALSE
  sample_freq <- 1
  set.seed(1)
  island_replicates <- list()
  out <- list()
  out[[1]] <- DAISIE:::DAISIE_sim_min_type2(
    time = totaltime,
    M = M,
    pars = pars,
    replicates = replicates,
    prop_type2_pool = prop_type2_pool,
    verbose = FALSE)
  island_replicates[[1]] <- out
  expect_output(
    formated_CS_sim <- DAISIE:::DAISIE_format_CS(
      island_replicates = island_replicates,
      time = time,
      M = mainland_n,
      sample_freq = sample_freq,
      island_type = island_type,
      verbose = verbose
    )
  )
})

test_that("abuse", {
  expect_error(DAISIE:::DAISIE_format_CS("nonsense"))
})
