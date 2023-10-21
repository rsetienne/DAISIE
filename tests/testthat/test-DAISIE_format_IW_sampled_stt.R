test_that("sampled stt, 1 type, no geodynamics, oceanic island (same arguments as geodynamics, 5 pars)", {
  pars <- c(0.5, 0.1, 10, 0.01, 0.5)
  total_time <- 1
  mainland_n <- 100
  verbose <- FALSE
  sample_freq <- 1
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
  replicates <- 1
  cond <- 1

  set.seed(1)
  island_replicates <- list()
  for (rep in seq_along(replicates)) {
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
  testthat::expect_silent(
    formatted_IW_sim <- DAISIE_format_IW_sampled_stt(
      island_replicates = island_replicates,
      total_time = total_time,
      M = mainland_n,
      sample_freq = sample_freq,
      verbose = verbose
    )
  )
  expected_IW_format <- list()
  expected_IW_format[[1]] <- list()
  stt_all <- matrix(ncol = 4, nrow = 2)
  colnames(stt_all) <- c("Time", "nI", "nA", "nC")
  stt_all[1, ] <- c(1, 0, 0, 0)
  stt_all[2, ] <- c(0, 0, 0, 3)
  expected_brts_table <- matrix(
    data = c(
      1.0, 0.0, 0.0, NA, NA,
      0.244818166871655, 1, 1, 1, NA,
      0.173128288990374, 1, 2, 1, NA,
      0.029668240213840, 1, 3, 1, NA
    ),
    ncol = 5,
    byrow = TRUE
  )
  colnames(expected_brts_table) <- c("brt", "clade", "event", "endemic", "col")

  expected_IW_format[[1]][[1]] <- list(island_age = 1,
                                       not_present = 99,
                                       stt_all = stt_all,
                                       brts_table = expected_brts_table)
  expected_IW_format[[1]][[2]] <- list(branching_times = c(1.000000000000000,
                                                           0.244818166871655,
                                                           0.173128288990374,
                                                           0.029668240213840),
                                       stac = 2,
                                       missing_species = 0)
  testthat::expect_equal(object = formatted_IW_sim, expected = expected_IW_format)

})

test_that("sampled stt, 1 type, geodynamics, oceanic island (same arguments as no geodynamics, 5 pars)", {
  total_time <- 1
  mainland_n <- 100
  verbose <- FALSE
  sample_freq <- 25
  set.seed(3)
  island_replicates <- list()
  out <- list()
  ext_pars <- c(1, 100)
  island_ontogeny <- 1
  sea_level <- 0
  pars <- c(0.0001, 2.2, 0.005, 1, 1)
  area_pars <- create_area_pars(
    max_area = 5000,
    current_area = 2500,
    proportional_peak_t = 0.5,
    total_island_age = 5,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0
  )
  hyper_pars <- create_hyper_pars(d = 0.2, x = 0.1)
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
  nonoceanic_pars <- c(0, 0)
  island_replicates[[1]] <- DAISIE_sim_core_time_dep(
    time = total_time,
    pars = pars,
    mainland_n = mainland_n,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    hyper_pars = hyper_pars,
    area_pars = area_pars,
    peak = peak,
    Amax = Amax,
    Amin = Amin,
    nonoceanic_pars = nonoceanic_pars
  )
  testthat::expect_silent(
    formatted_IW_sim <- DAISIE_format_IW_sampled_stt(
      island_replicates = island_replicates,
      total_time = total_time,
      M = mainland_n,
      sample_freq = sample_freq,
      verbose = verbose
    )
  )

  testthat::expect_equal(
    formatted_IW_sim[[1]][[1]]$island_age,
    total_time
  )
  testthat::expect_equal(
    formatted_IW_sim[[1]][[1]]$not_present,
    90
  )
  testthat::expect_equal(
    formatted_IW_sim[[1]][[1]]$stt_all[5, ],
    c(Time = 0.84, nI = 1.0, nA = 0.0, nC = 0.0)
  )
  testthat::expect_equal(
    formatted_IW_sim[[1]][[1]]$stt_all[12, ],
    c(Time = 0.56, nI = 4.0, nA = 1.0, nC = 0.0)
  )
  testthat::expect_equal(
    formatted_IW_sim[[1]][[1]]$stt_all[25, ],
    c(Time = 0.04, nI = 7.0, nA = 2.0, nC = 0.0)
  )

  testthat::expect_equal(
    formatted_IW_sim[[1]][[2]]$branching_times,
    c(1.00000000000000, 0.982693734281605)
  )

  testthat::expect_equal(
    formatted_IW_sim[[1]][[2]]$stac,
    2
  )

  testthat::expect_equal(
    formatted_IW_sim[[1]][[2]]$missing_species,
    0
  )
})

test_that("sampled stt, 1 type, no geodynamics, nonoceanic (same arguments as geodynamics, 5 pars)", {
  total_time <- 1
  clado_rate <- 2.550687345 # cladogenesis rate
  ext_rate <- 2.683454548 # extinction rate
  clade_carr_cap <- 10.0  # clade-level carrying capacity
  imm_rate <- 0.00933207 # immigration rate
  ana_rate <- 1.010073119 # anagenesis rate
  pars <- c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate)
  area_pars <- create_area_pars(
    max_area = 1,
    current_area = 1,
    proportional_peak_t = 0,
    total_island_age = 0,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0)
  hyper_pars <- create_hyper_pars(d = 0, x = 0)
  nonoceanic_pars <- c(0.1, 0.9)
  sample_freq <- 25
  mainland_n <- 100
  verbose <- FALSE
  set.seed(1)
  island_replicates <- list()
  island_replicates[[1]] <- DAISIE_sim_core_cr(
    time = total_time,
    pars = pars,
    mainland_n = mainland_n,
    area_pars = area_pars,
    hyper_pars = hyper_pars,
    nonoceanic_pars = nonoceanic_pars
  )
  testthat::expect_silent(
    formatted_IW_sim <- DAISIE_format_IW_sampled_stt(
      island_replicates = island_replicates,
      total_time = total_time,
      M = mainland_n,
      sample_freq = sample_freq,
      verbose = verbose
    )
  )
  testthat::expect_gt(sum(formatted_IW_sim[[1]][[1]]$stt_all[1, ]), total_time)
})
