context("DAISIE_format_CS")

test_that("silent with empty island with correct output", {
  pars <- c(0.4, 0.2, 10, 0.0001, 0.5)
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
    formatted_CS_sim <- DAISIE:::DAISIE_format_CS(
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
                                  missing_species = 0)
  expect_identical(formatted_CS_sim, expected_CS_format)
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
    formatted_CS_sim <- DAISIE:::DAISIE_format_CS(
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
                                       missing_species = 0)
  expect_equal(formatted_CS_sim, expected_CS_format)
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
    formatted_CS_sim <- DAISIE:::DAISIE_format_CS(
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
  skip("this scenario can't run (min_type2 always produces species)")
  pars <- c(0, 10, 1, 0.0001, 0, 0, 10, 1, 0.0001, 0)
  totaltime <- 1
  M <- 1
  mainland_n <- M
  replicates <- 1
  prop_type2_pool <- 0.1
  verbose <- FALSE
  sample_freq <- 1
  island_type <- "oceanic"
  set.seed(1)
  island_replicates <- list()
  island_replicates <- DAISIE:::DAISIE_sim_min_type2(
    time = totaltime,
    M = M,
    pars = pars,
    replicates = replicates,
    prop_type2_pool = prop_type2_pool,
    verbose = verbose
    )
  expect_silent(
    formatted_CS_sim <- DAISIE:::DAISIE_format_CS( # sim_min_type2 produces list with one extra element and fails
      island_replicates = island_replicates,
      time = totaltime,
      M = mainland_n,
      sample_freq = sample_freq,
      island_type = island_type,
      verbose = verbose
    )
  )
})

test_that("silent with non-empty 2 type island", {
  pars <- c(0.4, 0.1, 10, 1, 0.5, 0.4, 0.1, 10, 1, 0.5)
  totaltime <- 1
  M <- 10
  mainland_n <- M
  island_type <- "oceanic"
  verbose <- FALSE
  replicates <- 1
  sample_freq <- 10
  set.seed(1)
  island_replicates <- list()
  prop_type2_pool <- 0.4
  island_replicates <- DAISIE:::DAISIE_sim_min_type2(
    time = totaltime,
    M = M,
    pars = pars,
    replicates = replicates,
    prop_type2_pool = prop_type2_pool,
    verbose = FALSE)
  expect_silent(
    formatted_CS_sim <- DAISIE:::DAISIE_format_CS(
      island_replicates = island_replicates,
      time = totaltime,
      M = mainland_n,
      sample_freq = sample_freq,
      island_type = island_type,
      verbose = verbose
    )
  )
})

test_that("silent with non-empty 2 type island full stt", {
  pars <- c(0.4, 0.1, 10, 1, 0.5, 0.4, 0.1, 10, 1, 0.5)
  totaltime <- 1
  M <- 10
  mainland_n <- M
  island_type <- "oceanic"
  verbose <- FALSE
  replicates <- 1
  sample_freq <- Inf
  set.seed(1)
  island_replicates <- list()
  prop_type2_pool <- 0.4
  island_replicates <- DAISIE:::DAISIE_sim_min_type2(
    time = totaltime,
    M = M,
    pars = pars,
    replicates = replicates,
    prop_type2_pool = prop_type2_pool,
    verbose = FALSE)
  expect_silent(
    formatted_CS_sim <- DAISIE:::DAISIE_format_CS(
      island_replicates = island_replicates,
      time = totaltime,
      M = mainland_n,
      island_type = island_type,
      verbose = verbose
    )
  )
})

test_that("abuse", {
  expect_error(
    DAISIE:::DAISIE_format_CS(
      island_replicates = "nonsense",
      time = "nonsense",
      M = "nonsense",
      sample_freq = "nonsense",
      island_type = "nonsense",
      verbose = "nonsense"
    )
  )
})


test_that("use full stt", {
  pars <- c(0.4, 0.2, 10, 2, 0.5)
  time <- 5
  mainland_n <- 1
  verbose <- FALSE
  sample_freq <- Inf
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
    formatted_CS_sim <- DAISIE:::DAISIE_format_CS(
      island_replicates = island_replicates,
      time = time,
      M = mainland_n,
      sample_freq = sample_freq,
      island_type = island_type,
      verbose = verbose
    )
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$island_age,
    5
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$not_present,
    0
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$stt_all[12, ],
    c(Time = 1.5716508023714537, nI = 0.0, nA = 0.0, nC = 2.0, present = 2.0)
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$stt_all[5, ],
    c(Time = 3.8154825768724887, nI = 0.0, nA = 1.0, nC = 0.0, present = 1.0)
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$stt_all[19, ],
    c(Time = 0.09210138119067679, nI = 1.0, nA = 1.0, nC = 2.0, present = 4.0)
  )

  expect_equal(
    formatted_CS_sim[[1]][[2]]$branching_times,
    c(5.0, 1.348741816972570007, 0.092101381190680301)
  )

  expect_equal(
    formatted_CS_sim[[1]][[2]]$stac,
    3
  )

  expect_equal(
    formatted_CS_sim[[1]][[2]]$missing_species,
    0
  )
})

test_that("use complete stt with ontogeny", {
  totaltime <- 10
  mainland_n <- 1
  verbose <- FALSE
  sample_freq <- Inf
  set.seed(1)
  island_replicates <- list()
  out <- list()

  pars = c(0.0001, 2.2, 0.005, 1, 1)
  island_type = "oceanic"
  area_pars = create_area_pars(
    max_area = 5000,
    proportional_peak_t = 0.5,
    peak_sharpness = 1,
    total_island_age = 15,
    sea_level_amplitude = 0,
    sea_level_frequency = 0
  )
  ext_pars = c(1, 100)
  island_ontogeny = 1
  sea_level = "const"
  out[[1]] <- DAISIE:::DAISIE_sim_core(
    time = totaltime,
    pars = pars,
    mainland_n = mainland_n,
    island_ontogeny = island_ontogeny,
    area_pars = area_pars,
    ext_pars = ext_pars,
    island_type = island_type,
    sea_level = sea_level
  )
  island_replicates[[1]] <- out
  expect_silent(
    formatted_CS_sim <- DAISIE:::DAISIE_format_CS(
      island_replicates = island_replicates,
      time = totaltime,
      M = mainland_n,
      sample_freq = sample_freq,
      island_type = island_type,
      verbose = verbose
    )
  )

  expect_equal(
    formatted_CS_sim[[1]][[1]]$island_age,
    10
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$not_present,
    0
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$stt_all[12, ],
    c(Time = 5.6629724151660916, nI = 1.0, nA = 1.0, nC = 0.0, present = 2.0)
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$stt_all[5, ],
    c(Time = 7.3934919638882635, nI = 1.0, nA = 0.0, nC = 0.0, present = 1.0)
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$stt_all[25, ],
    c(Time = 0.85034199874260885, nI = 0.0, nA = 1.0, nC = 0.0, present = 1.0)
  )

  expect_equal(
    formatted_CS_sim[[1]][[2]]$branching_times,
    c(10.0, 0.67565477313507005)
  )

  expect_equal(
    formatted_CS_sim[[1]][[2]]$stac,
    2
  )

  expect_equal(
    formatted_CS_sim[[1]][[2]]$missing_species,
    0
  )
})

test_that("full stt works with multiple replicates", {
  pars <- c(0.4, 0.2, 10, 2, 0.5)
  time <- 5
  mainland_n <- 1
  verbose <- FALSE
  sample_freq <- Inf
  island_type <- "oceanic"
  set.seed(1)
  island_replicates <- list()
  island_replicates[[1]] <- island_replicates
  out <- list()
  out[[1]] <- DAISIE:::DAISIE_sim_core(
    time = time,
    pars = pars,
    mainland_n = mainland_n
  )
  out[[2]] <- DAISIE:::DAISIE_sim_core(
    time = time,
    pars = pars,
    mainland_n = mainland_n
  )
  island_replicates[[1]] <- out
  expect_silent(
    formatted_CS_sim <- DAISIE:::DAISIE_format_CS(
      island_replicates = island_replicates,
      time = time,
      M = mainland_n,
      sample_freq = sample_freq,
      island_type = island_type,
      verbose = verbose
    )
  )
})

test_that("full stt works with empty island", {
  pars <- c(0.4, 0.2, 10, 0.0000001, 0.5)
  totaltime <- 1
  mainland_n <- 1
  verbose <- FALSE
  sample_freq <- Inf
  island_type <- "oceanic"
  set.seed(1)
  replicates <- 2
  island_replicates <- list()
  island_ontogeny <- 0
  for (rep in 1:replicates) {
    island_replicates[[rep]] <- list()
    full_list <- list()
    out <- list()
    for (m_spec in 1:mainland_n) {
      out$branching_times <- c(10)
        out <- DAISIE:::DAISIE_sim_core(
          island_ontogeny = island_ontogeny,
          time = totaltime,
          mainland_n = 1,
          pars = pars,
          island_type = island_type
        )
      full_list[[m_spec]] <- out
    }
    island_replicates[[rep]] <- full_list
  }
  expect_silent(
    formatted_CS_sim <- DAISIE:::DAISIE_format_CS(
      island_replicates = island_replicates,
      time = time,
      M = mainland_n,
      sample_freq = sample_freq,
      island_type = island_type,
      verbose = verbose
    )
  )
})
