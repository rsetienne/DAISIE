context("DAISIE_format_CS_full_stt")

test_that("complete stt, 1 type, no geodynamics, oceanic island (same arguments
          as geodynamics, 5 pars)", {
  pars <- c(0.4, 0.2, 10, 2, 0.5)
  totaltime <- 1
  mainland_n <- 2
  verbose <- FALSE
  set.seed(1)
  replicates <- 3
  island_replicates <- list()
  for (rep in 1:replicates) {
    island_replicates[[rep]] <- list()
    full_list <- list()
    out <- list()
    for (m_spec in 1:mainland_n) {
      out$branching_times <- c(10)
      while (length(out$branching_times) == 1) {
        out <- DAISIE:::DAISIE_sim_core_constant_rate(
          time = totaltime,
          mainland_n = 1,
          pars = pars
        )
      }
      full_list[[m_spec]] <- out
    }
    island_replicates[[rep]] <- full_list
  }


  expect_silent(
    formatted_CS_sim <- DAISIE:::DAISIE_format_CS_full_stt(
      island_replicates = island_replicates,
      time = totaltime,
      M = mainland_n,
      verbose = verbose
    )
  )

  expect_equal(
    formatted_CS_sim[[1]][[1]]$island_age,
    1
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$not_present,
    0
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$stt_all[2, ],
    c(Time = 0.8129771533422172, nI = 1.0, nA = 0.0, nC = 0.0, present = 1.0)
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$stt_all[4, ],
    c(Time = 0.5735296212440449, nI = 1.0, nA = 0.0, nC = 0.0, present = 1.0)
  )

  expect_equal(
    formatted_CS_sim[[1]][[2]]$branching_times,
    c(1.000000000000000, 0.444415577345164)
  )

  expect_equal(
    formatted_CS_sim[[1]][[2]]$stac,
    4
  )

  expect_equal(
    formatted_CS_sim[[1]][[2]]$missing_species,
    0
  )

  expect_equal(
    formatted_CS_sim[[1]][[3]]$branching_times,
    c(1.000000000000000, 0.326839228620516)
  )

  expect_equal(
    formatted_CS_sim[[1]][[3]]$stac,
    4
  )

  expect_equal(
    formatted_CS_sim[[1]][[3]]$missing_species,
    0
  )
})

test_that("complete stt, 1 type, geodynamics, oceanic island (same arguments as
          no geodynamics, 5 pars)", {
  totaltime <- 5
  mainland_n <- 2
  verbose <- FALSE
  set.seed(1)
  island_replicates <- list()
  replicates <- 3
  island_ontogeny <- 1
  sea_level <- 0
  pars <- c(0.0001, 2.2, 0.005, 1, 1)
  default_pars <- create_default_pars(
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    area_pars = create_area_pars(
      max_area = 5000,
      proportional_peak_t = 0.5,
      peak_sharpness = 1,
      total_island_age = 15,
      sea_level_amplitude = 0,
      sea_level_frequency = 0,
      island_gradient_angle = 0
    ),
    ext_pars = c(1, 100),
    hyper_pars = NULL,
    dist_pars = NULL,
    totaltime = totaltime,
    pars = pars
  )

  for (rep in 1:replicates) {
    island_replicates[[rep]] <- list()
    full_list <- list()
    out <- list()
    for (m_spec in 1:mainland_n) {
      out$branching_times <- c(10)
      while (length(out$branching_times) == 1) {
        out <- DAISIE:::DAISIE_sim_core_time_dependent(
          island_ontogeny = 1,
          time = totaltime,
          mainland_n = 1,
          pars = pars,
          sea_level = sea_level,
          area_pars = default_pars$area_pars,
          ext_pars = default_pars$ext_pars,
          hyper_pars = default_pars$hyper_pars,
          extcutoff = 100,
          dist_pars = default_pars$dist_pars
        )
      }
      full_list[[m_spec]] <- out
    }
    island_replicates[[rep]] <- full_list
  }

  expect_silent(
    formatted_CS_sim <- DAISIE:::DAISIE_format_CS_full_stt(
      island_replicates = island_replicates,
      time = totaltime,
      M = mainland_n,
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
    formatted_CS_sim[[1]][[1]]$stt_all[2, ],
    c(Time = 4.4979918331304169, nI = 1.0, nA = 0.0, nC = 0.0, present = 1.0)
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$stt_all[5, ],
    c(Time = 4.2448181668716547, nI = 2.0, nA = 0.0, nC = 0.0, present = 2.0)
  )
  expect_equal(
    formatted_CS_sim[[1]][[2]]$branching_times,
    c(5.00000000000000, 0.69194184070825004, 0.53032188808254999)
  )
  expect_equal(
    formatted_CS_sim[[1]][[2]]$stac,
    2
  )
  expect_equal(
    formatted_CS_sim[[1]][[2]]$missing_species,
    0
  )

  # Correct number of rows on initial and final matrices
  expect_equal(
    nrow(formatted_CS_sim[[1]][[1]]$stt_all),
    (nrow(island_replicates[[1]][[1]]$stt_table) +
       nrow(island_replicates[[1]][[2]]$stt_table)) - 2
  )
  expect_equal(
    nrow(formatted_CS_sim[[2]][[1]]$stt_all),
    (nrow(island_replicates[[2]][[1]]$stt_table) +
       nrow(island_replicates[[2]][[2]]$stt_table)) - 2
  )
})

test_that("complete stt, 2 type, no geodynamics, oceanic island (same arguments
          as geodynamics, 10 pars)", {
  pars <- c(0.4, 0.1, 10, 1, 0.5, 0.4, 0.1, 10, 1, 0.5)
  totaltime <- 5
  M <- 10
  mainland_n <- M
  verbose <- FALSE
  replicates <- 2
  set.seed(1)
  island_replicates <- list()
  prop_type2_pool <- 0.4

  island_replicates <- DAISIE:::DAISIE_sim_min_type2(
    time = totaltime,
    M = M,
    pars = pars,
    replicates = replicates,
    prop_type2_pool = prop_type2_pool,
    verbose = FALSE
  )
  expect_silent(
    formatted_CS_sim <- DAISIE:::DAISIE_format_CS_full_stt(
      island_replicates = island_replicates,
      time = totaltime,
      M = mainland_n,
      verbose = verbose
    )
  )

  expect_equal(
    formatted_CS_sim[[1]][[1]]$island_age,
    5
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$not_present_type1,
    0
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$not_present_type2,
    0
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$stt_all[2, ],
    c(Time = 4.979749586898834, nI = 1.0, nA = 0.0, nC = 0.0, present = 1.0)
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$stt_all[5, ],
    c(Time = 4.590700004177261, nI = 1.0, nA = 0.0, nC = 2.0, present = 3.0)
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$stt_all[19, ],
    c(Time = 3.867811945494177, nI = 6.0, nA = 3.0, nC = 4.0, present = 13.0)
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$stt_type1[5, ],
    c(Time = 4.460097487432346, nI = 2.0, nA = 0.0, nC = 2.0, present = 4.0)
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$stt_type1[15, ],
    c(Time = 3.771978362895029, nI = 4.0, nA = 1.0, nC = 6.0, present = 11.0)
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$stt_type2[5, ],
    c(Time = 4.317296566441655, nI = 2.0, nA = 0.0, nC = 0.0, present = 2.0)
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$stt_type2[15, ],
    c(Time = 3.164011375762181, nI = 3.0, nA = 2.0, nC = 0.0, present = 5.0)
  )

  expect_equal(
    formatted_CS_sim[[1]][[2]]$branching_times,
    c(5.00000000000000, 4.24481816687165, 1.44735043895909)
  )

  expect_equal(
    formatted_CS_sim[[1]][[2]]$stac,
    3
  )

  expect_equal(
    formatted_CS_sim[[1]][[2]]$missing_species,
    0
  )

  # Correct number of rows on initial and final matrices. Make sure length
  # of full stt is the sum of individual colonist stt table,
  # minus the first and last lines of each stt table. The first and last lines
  # are always the same in all stt tables.
  stt_nrows <- sapply(island_replicates[[1]], FUN = function(x) nrow(x$stt_table))
  sum_stt_nrows <- sum(stt_nrows)
  expect_equal(
    nrow(formatted_CS_sim[[1]][[1]]$stt_all),
    sum_stt_nrows - 18
  )
  stt_nrows <- sapply(island_replicates[[2]], FUN = function(x) nrow(x$stt_table))
  sum_stt_nrows <- sum(stt_nrows)
  expect_equal(
    nrow(formatted_CS_sim[[2]][[1]]$stt_all),
    sum_stt_nrows - 18
  )

  # Type 1
  type_1s <- island_replicates[[1]][1:6]
  stt_nrows <- sapply(type_1s, FUN = function(x) nrow(x$stt_table))
  sum_stt_nrows <- sum(stt_nrows)
  expect_equal(
    nrow(formatted_CS_sim[[1]][[1]]$stt_type1),
    sum_stt_nrows - 10
  )
  type_1s <- island_replicates[[2]][1:6]
  stt_nrows <- sapply(type_1s, FUN = function(x) nrow(x$stt_table))
  sum_stt_nrows <- sum(stt_nrows)
  expect_equal(
    nrow(formatted_CS_sim[[2]][[1]]$stt_type1),
    sum_stt_nrows - 10
  )

  # Type 2
  type_2s <- island_replicates[[1]][7:10]
  stt_nrows <- sapply(type_2s, FUN = function(x) nrow(x$stt_table))
  sum_stt_nrows <- sum(stt_nrows)
  expect_equal(
    nrow(formatted_CS_sim[[1]][[1]]$stt_type2),
    sum_stt_nrows - 6
  )
  type_2s <- island_replicates[[2]][7:10]
  stt_nrows <- sapply(type_2s, FUN = function(x) nrow(x$stt_table))
  sum_stt_nrows <- sum(stt_nrows)
  expect_equal(
    nrow(formatted_CS_sim[[2]][[1]]$stt_type2),
    sum_stt_nrows - 6
  )
})

test_that("complete stt, 1 type, no geodynamics, nonoceanic (same arguments as geodynamics, 5 pars)", {
  totaltime <- 3
  mainland_n <- 2
  clado_rate <- 1 # cladogenesis rate
  ext_rate <- 0.3 # extinction rate
  clade_carr_cap <- 10.0  # clade-level carrying capacity
  imm_rate <- 0.00933207 # immigration rate
  ana_rate <- 1.010073119 # anagenesis rate
  pars <- c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate)
  replicates <- 3
  nonoceanic_pars <- c(0.1, 0.9)
  island_replicates <- list()
  verbose <- FALSE

  for (rep in 1:replicates) {
    island_replicates[[rep]] <- list()
    full_list <- list()
    out <- list()
    for (m_spec in 1:mainland_n) {
      out$branching_times <- c(10)
      while (length(out$branching_times) == 1) {
        out <- DAISIE:::DAISIE_sim_core_constant_rate(
          time = totaltime,
          mainland_n = 1,
          pars = pars,
          nonoceanic_pars = nonoceanic_pars
        )
      }
      full_list[[m_spec]] <- out
    }
    island_replicates[[rep]] <- full_list
  }
  expect_silent(
    formatted_CS_sim <- DAISIE:::DAISIE_format_CS_full_stt(
      island_replicates = island_replicates,
      time = totaltime,
      M = mainland_n,
      verbose = verbose
    )
  )
})

test_that("complete stt, 1 type, no geodynamics, oceanic island (same
          arguments as geodynamics, 5 pars) verbose", {
  pars <- c(0.4, 0.2, 10, 2, 0.8)
  totaltime <- 1
  mainland_n <- 2
  verbose <- TRUE
  set.seed(1)
  replicates <- 3
  island_ontogeny <- 0

  island_replicates <- list()
  for (rep in 1:replicates) {
    island_replicates[[rep]] <- list()
    full_list <- list()
    out <- list()
    for (m_spec in 1:mainland_n) {
      out$branching_times <- c(10)
      while (length(out$branching_times) == 1) {
        out <- DAISIE:::DAISIE_sim_core_constant_rate(
          time = totaltime,
          mainland_n = 1,
          pars = pars
        )
      }
      full_list[[m_spec]] <- out
    }
    island_replicates[[rep]] <- full_list
  }
  expect_output(
    formatted_CS_sim <- DAISIE:::DAISIE_format_CS_full_stt(
      island_replicates = island_replicates,
      time = totaltime,
      M = mainland_n,
      verbose = verbose
    ),
    regexp = "Island being formatted: 3/3"
  )
})

