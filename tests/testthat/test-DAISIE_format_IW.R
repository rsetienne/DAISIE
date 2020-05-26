context("test-DAISIE_format_IW")

test_that("use with empty island", {
  pars <- c(0.4, 0.2, 10, 0.0001, 0.5)
  time <- 1
  mainland_n <- 1000
  verbose <- FALSE
  sample_freq <- 1
  start_midway <- FALSE
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
      verbose = verbose
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
      verbose = verbose
    )
  )
})

test_that("abuse", {
  expect_error(DAISIE:::DAISIE_format_IW(
    island_replicates = "nonsense",
    time = time,
    M = mainland_n,
    sample_freq = sample_freq,
    verbose = verbose
  ))
})

test_that("new and v1.5 give same results", {

  tol <- 1e-13
  sim_time <- 10
  n_mainland_species <- 300
  clado_rate <- 1.0
  ext_rate <- 0.5
  carr_cap <- 10
  imm_rate <- 1.0
  ana_rate <- 1.0
  pars <- c(clado_rate, ext_rate, carr_cap, imm_rate, ana_rate)
  rng_seed <- 42
  set.seed(rng_seed)
  ff <- DAISIE_sim(time = sim_time,
                   M = n_mainland_species,
                   pars = pars,
                   replicates = 1,
                   divdepmodel = 'IW')
  ff[[1]][[1]]$brts_table <- 0
  new <- DAISIE:::Add_brt_table(ff[[1]])
  new <- new[[1]]$brts_table[-1,]
  old <- DAISIE:::Add_brt_table_v1_5(ff[[1]])
  old <- old[[1]]$brts_table[-1,]
  testthat::expect_true(all(abs(new - old) < tol))
})
