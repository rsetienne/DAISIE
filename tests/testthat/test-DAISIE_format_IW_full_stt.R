test_that("complete stt, 1 type, no geodynamics, oceanic island, one trait state
          (same arguments as geodynamics, 5 pars)", {
  pars <- c(0.6, 0.1, 10, 0.01, 0.5)
  total_time <- 1
  mainland_n <- 100
  area_pars <- create_area_pars(
    max_area = 1,
    current_area = 1,
    proportional_peak_t = 0,
    total_island_age = 0,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0)
  hyper_pars <- create_hyper_pars(d = 0, x = 0)
  nonoceanic_pars <- c(0, 0)
  verbose <- FALSE
  set.seed(1)
  replicates <- 2
  cond <- 1
  island_replicates <- list()

  for (rep in seq_len(replicates)) {
    if (cond == 0) {
      number_present <- -1
    } else {
      number_present <- 0
    }
    while (number_present < cond) {
      island_replicates[[rep]] <- DAISIE_sim_core_cr(
        time = total_time,
        mainland_n = mainland_n,
        pars = pars,
        nonoceanic_pars = nonoceanic_pars,
        hyper_pars = hyper_pars,
        area_pars = area_pars
      )
      stac_vec <- unlist(island_replicates)[which(
        names(unlist(island_replicates)) == "taxon_list.stac"
      )]
      present <- which(stac_vec != 0)
      number_present <- length(present)
    }
  }
  expect_silent(
    formatted_IW_sim <- DAISIE_format_IW_full_stt(
      island_replicates = island_replicates,
      total_time = total_time,
      M = mainland_n,
      verbose = verbose
    )
  )

  expect_equal(
    formatted_IW_sim[[1]][[1]]$island_age,
    1
  )
  expect_equal(
    formatted_IW_sim[[1]][[1]]$not_present,
    99
  )
  expect_equal(
    formatted_IW_sim[[1]][[1]]$stt_all[2, ],
    c(Time = 0.244818166871655, nI = 1.0, nA = 0.0, nC = 0.0)
  )
  expect_equal(
    formatted_IW_sim[[1]][[1]]$stt_all[4, ],
    c(Time = 0.0, nI = 0.0, nA = 1.0, nC = 0.0)
  )
  expected_brts_table <- matrix(
    c(1.0, 0.24481816687165, 0, 1, 0, 1, NA, 1, NA, NA),
    nrow = 2
  )
  colnames(expected_brts_table) <- c("brt", "clade", "event", "endemic", "col")
  rownames(expected_brts_table) <- c("", "brts_table")
  expect_equal(
    formatted_IW_sim[[1]][[1]]$brts_table, expected_brts_table
  )

  expect_equal(
    formatted_IW_sim[[1]][[2]]$branching_times,
    c(1.000000000000000, 0.244818166871655)
  )

  expect_equal(
    formatted_IW_sim[[1]][[2]]$stac,
    2
  )

  expect_equal(
    formatted_IW_sim[[1]][[2]]$missing_species,
    0
  )

  expect_equal(
    formatted_IW_sim[[2]][[1]]$island_age,
    1
  )
  expect_equal(
    formatted_IW_sim[[2]][[1]]$not_present,
    99
  )
  expect_equal(
    formatted_IW_sim[[2]][[1]]$stt_all[2, ],
    c(Time = 0.741771912202239, nI = 1.0, nA = 0.0, nC = 0.0)
  )
  expect_equal(
    formatted_IW_sim[[2]][[1]]$stt_all[3, ],
    c(Time = 0.0, nI = 1.0, nA = 0.0, nC = 0.0)
  )

  expected_brts_table <- matrix(
    c(1.0, 0.741771912202239, 0, 1, 0, 1, NA, 0, NA, NA),
    nrow = 2
  )
  colnames(expected_brts_table) <- c("brt", "clade", "event", "endemic", "col")
  rownames(expected_brts_table) <- c("", "brts_table")
  expect_equal(
    formatted_IW_sim[[2]][[1]]$brts_table, expected_brts_table
  )


  expect_equal(
    formatted_IW_sim[[2]][[2]]$branching_times,
    c(1.000000000000000, 0.741771912202239)
  )

  expect_equal(
    formatted_IW_sim[[2]][[2]]$stac,
    4
  )

  expect_equal(
    formatted_IW_sim[[2]][[2]]$missing_species,
    0
  )
})


test_that("complete stt, 1 type, geodynamics, oceanic island, one trait state
          (same arguments as no geodynamics, 5 pars)", {
  total_time <- 5
  mainland_n <- 100
  verbose <- FALSE
  set.seed(1)
  replicates <- 2
  island_ontogeny <- 1
  sea_level <- 0
  pars <- c(0.0001, 2.2, 0.005, 1, 1)
  area_pars <- create_area_pars(
    max_area = 5000,
    current_area = 2500,
    proportional_peak_t = 0.5,
    total_island_age = 15,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0
  )
  hyper_pars <- create_hyper_pars(d = 0.2, x = 0.1)
  nonoceanic_pars <- c(0, 0)
  peak <- calc_peak(total_time = total_time,
                    area_pars = area_pars)
  Amax <- get_global_max_area(total_time = total_time,
                              area_pars = area_pars,
                              peak = peak,
                              island_ontogeny = island_ontogeny,
                              sea_level = sea_level)
  Amin <- get_global_min_area(total_time = total_time,
                              area_pars = area_pars,
                              peak = peak,
                              island_ontogeny = island_ontogeny,
                              sea_level = sea_level)
  cond <- 1

  island_replicates <- list()

  for (rep in seq_len(replicates)) {
    if (cond == 0) {
      number_present <- -1
    } else {
      number_present <- 0
    }
    while (number_present < cond) {
      island_replicates[[rep]] <- DAISIE_sim_core_time_dep(
        time = total_time,
        mainland_n = mainland_n,
        pars = pars,
        nonoceanic_pars = nonoceanic_pars,
        hyper_pars = hyper_pars,
        area_pars = area_pars,
        peak = peak,
        Amax = Amax,
        Amin = Amin,
        island_ontogeny = island_ontogeny,
        sea_level = sea_level
      )
      stac_vec <- unlist(island_replicates)[which(
        names(unlist(island_replicates)) == "taxon_list.stac"
      )]
      present <- which(stac_vec != 0)
      number_present <- length(present)
    }
  }

  expect_silent(
    formatted_IW_sim <- DAISIE_format_IW_full_stt(
      island_replicates = island_replicates,
      total_time = total_time,
      M = mainland_n,
      verbose = verbose
    )
  )

  expect_equal(
    formatted_IW_sim[[1]][[1]]$island_age,
    5
  )
  expect_equal(
    formatted_IW_sim[[1]][[1]]$not_present,
    90
  )
  expect_equal(
    formatted_IW_sim[[1]][[1]]$stt_all[2, ],
    c(Time = 4.99244818166872, nI = 1.0, nA = 0.0, nC = 0.0)
  )
  expect_equal(
    formatted_IW_sim[[1]][[1]]$stt_all[5, ],
    c(Time = 4.79587602350116, nI = 0.0, nA = 0.0, nC = 0.0)
  )
  expect_equal(
    formatted_IW_sim[[1]][[2]]$branching_times,
    c(5.00000000000000, 1.38252473752213)
  )
  expect_equal(
    formatted_IW_sim[[1]][[2]]$stac,
    4
  )
  expect_equal(
    formatted_IW_sim[[1]][[2]]$missing_species,
    0
  )

  expect_equal(
    formatted_IW_sim[[1]][[1]]$brts_table[1, ],
    c(brt = 5, clade = 0, event = 0, endemic = NA, col = NA)
  )
  expect_equal(
    formatted_IW_sim[[1]][[1]]$brts_table[5, ],
    c(brt = 0.83094531417507, clade = 4, event = 1, endemic = 1, col = NA)
  )

})


test_that("complete stt, 1 type, no geodynamics, nonoceanic, one trait state
          (same arguments as geodynamics, 5 pars)", {
  total_time <- 3
  mainland_n <- 100
  clado_rate <- 1 # cladogenesis rate
  ext_rate <- 0.3 # extinction rate
  clade_carr_cap <- 10.0  # clade-level carrying capacity
  imm_rate <- 0.00933207 # immigration rate
  ana_rate <- 1.010073119 # anagenesis rate
  pars <- c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate)
  replicates <- 2
  hyper_pars <- create_hyper_pars(d = 0, x = 0)
  nonoceanic_pars <- c(0.1, 0.9)
  verbose <- FALSE
  area_pars <- create_area_pars(
    max_area = 1,
    current_area = 1,
    proportional_peak_t = 0,
    total_island_age = 0,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0)
  hyper_pars <- create_hyper_pars(d = 0, x = 0)

  island_replicates <- list()

  for (rep in seq_len(replicates)) {
    if (cond == 0) {
      number_present <- -1
    } else {
      number_present <- 0
    }
    while (number_present < cond) {
      island_replicates[[rep]] <- DAISIE_sim_core_cr(
        time = total_time,
        mainland_n = mainland_n,
        pars = pars,
        area_pars = area_pars,
        hyper_pars = hyper_pars,
        nonoceanic_pars = nonoceanic_pars
      )
      stac_vec <- unlist(island_replicates)[which(
        names(unlist(island_replicates)) == "taxon_list.stac"
      )]
      present <- which(stac_vec != 0)
      number_present <- length(present)
    }
  }
  # stt is not empty at start
  expect_gte(sum(island_replicates[[1]]$stt_table[1, ]), expected = total_time)
  expect_gte(sum(island_replicates[[2]]$stt_table[1, ]), expected = total_time)
  expect_silent(
    formatted_IW_sim <- DAISIE_format_IW_full_stt(
      island_replicates = island_replicates,
      total_time = total_time,
      M = mainland_n,
      verbose = verbose
    )
  )
})

test_that("complete stt, 1 type, no geodynamics, oceanic island, one trait
          state (same arguments as geodynamics, 5 pars) verbose", {
  pars <- c(0.4, 0.2, 10, 2, 0.8)
  total_time <- 1
  mainland_n <- 100
  verbose <- TRUE
  area_pars <- create_area_pars(
    max_area = 1,
    current_area = 1,
    proportional_peak_t = 0,
    total_island_age = 0,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0)
  hyper_pars <- create_hyper_pars(d = 0, x = 0)
  nonoceanic_pars <- c(0, 0)
  set.seed(1)
  replicates <- 2
  island_ontogeny <- 0
  sea_level <- 0
  cond <- 1

  for (rep in seq_len(replicates)) {
    if (cond == 0) {
      number_present <- -1
    } else {
      number_present <- 0
    }
    while (number_present < cond) {
      island_replicates[[rep]] <- DAISIE_sim_core_cr(
        time = total_time,
        mainland_n = mainland_n,
        pars = pars,
        nonoceanic_pars = nonoceanic_pars,
        hyper_pars = hyper_pars,
        area_pars = area_pars
      )
      stac_vec <- unlist(island_replicates)[which(
        names(unlist(island_replicates)) == "taxon_list.stac"
      )]
      present <- which(stac_vec != 0)
      number_present <- length(present)
    }
  }
  expect_message(
    formatted_IW_sim <- DAISIE_format_IW_full_stt(
      island_replicates = island_replicates,
      total_time = total_time,
      M = mainland_n,
      verbose = verbose
    ),
    regexp = "Island being formatted:*"
  )
})


test_that("when no colonization happens returns 0", {
  pars <- c(0.4, 0.2, 10, 0.000001, 0.5)
  total_time <- 1
  mainland_n <- 100
  area_pars <- create_area_pars(
    max_area = 1,
    current_area = 1,
    proportional_peak_t = 0,
    total_island_age = 0,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0)
  hyper_pars <- create_hyper_pars(d = 0, x = 0)
  nonoceanic_pars <- c(0, 0)
  verbose <- FALSE
  set.seed(2)
  replicates <- 1
  island_replicates <- list()
  island_replicates[[1]] <- list()
  out <- DAISIE_sim_core_cr(
    time = total_time,
    mainland_n = 1,
    pars = pars,
    area_pars = area_pars,
    hyper_pars = hyper_pars,
    nonoceanic_pars = nonoceanic_pars
  )


  island_replicates <- list(out)

  expect_silent(
    formatted_IW_sim <- DAISIE_format_IW_full_stt(
      island_replicates = island_replicates,
      total_time = total_time,
      M = mainland_n,
      verbose = verbose
    )
  )

  expect_equal(
    formatted_IW_sim[[1]][[1]]$island_age,
    1
  )
  expect_equal(
    formatted_IW_sim[[1]][[1]]$not_present,
    100
  )
  expect_equal(
    formatted_IW_sim[[1]][[1]]$stt_all[2, ],
    c(Time = 0, nI = 0.0, nA = 0.0, nC = 0.0)
  )
})




