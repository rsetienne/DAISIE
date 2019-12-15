context("DAISIE_format_CS")

test_that("sampled stt, 1 type, no geodynamics, oceanic island (same arguments as geodynamics, 5 pars)", {
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
                                       missing_species = 0)
  expect_equal(formated_CS_sim, expected_CS_format)
})
test_that("sampled stt, 1 type, geodynamics, oceanic island (same arguments as no geodynamics, 5 pars)", {
  time <- 10
  mainland_n <- 1
  verbose <- FALSE
  sample_freq <- 25
  set.seed(1)
  island_replicates <- list()
  out <- list()

  pars = c(0.0001, 2.2, 0.005, 1, 1)
  ddmodel_sim = 11
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
    time = time,
    pars = pars,
    mainland_n = mainland_n,
    island_ontogeny = island_ontogeny,
    area_pars = area_pars,
    ext_pars = ext_pars,
    ddmodel_sim = ddmodel_sim,
    island_type = island_type,
    sea_level = sea_level
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

  expect_equal(
    formated_CS_sim[[1]][[1]]$island_age,
    10
  )
  expect_equal(
    formated_CS_sim[[1]][[1]]$not_present,
    0
  )
  expect_equal(
    formated_CS_sim[[1]][[1]]$stt_all[12, ],
    c(Time = 5.600000000000001, nI = 1.0, nA = 1.0, nC = 0.0, present = 1.0)
  )
  expect_equal(
    formated_CS_sim[[1]][[1]]$stt_all[5, ],
    c(Time = 8.4, nI = 1.0, nA = 0.0, nC = 0.0, present = 1.0)
  )
  expect_equal(
    formated_CS_sim[[1]][[1]]$stt_all[25, ],
    c(Time = 0.4, nI = 0.0, nA = 1.0, nC = 0.0, present = 1.0)
  )

  expect_equal(
    formated_CS_sim[[1]][[2]]$branching_times,
    c(10.0, 0.67565477313507005)
  )

  expect_equal(
    formated_CS_sim[[1]][[2]]$stac,
    2
  )

  expect_equal(
    formated_CS_sim[[1]][[2]]$missing_species,
    0
  )
})

test_that("sampled stt, 2 type, no geodynamics, oceanic island (same arguments as geodynamics, 10 pars)", {
  pars <- c(0.4, 0.1, 10, 1, 0.5, 0.4, 0.1, 10, 1, 0.5)
  totaltime <- 1
  M <- 10
  mainland_n <- M
  island_type <- "oceanic"
  verbose <- FALSE
  replicates <- 1
  sample_freq <- 25
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
    formated_CS_sim <- DAISIE:::DAISIE_format_CS(
      island_replicates = island_replicates,
      time = totaltime,
      M = mainland_n,
      sample_freq = sample_freq,
      island_type = island_type,
      verbose = verbose
    )
  )

  expect_equal(
    names(formated_CS_sim[[1]][[1]]),
    c(
      "island_age",
      "not_present_type1",
      "not_present_type2",
      "stt_all",
      "stt_type1",
      "stt_type2"
    )
  )

  expect_equal(
    names(formated_CS_sim[[1]][[2]]),
    c(
      "branching_times",
      "stac",
      "missing_species",
      "other_clades_same_ancestor",
      "type1or2"
    )
  )

  expect_equal(
    length(formated_CS_sim[[1]]), 8
  )

  # Sampled STT has the correct size
  expect_equal(
    dim(formated_CS_sim[[1]][[1]]$stt_all), c(26, 5)
  )
  expect_equal(
    dim(formated_CS_sim[[1]][[1]]$stt_type1), c(26, 5)
  )
  expect_equal(
    dim(formated_CS_sim[[1]][[1]]$stt_type2), c(26, 5)
  )

})
test_that("sampled stt, 2 type, geodynamics, oceanic island (same arguments as geodynamics, 10 pars)", {
  skip("DAISIE_sim_min_type2 can't run with geodynamics")
  pars <- c(0.4, 0.1, 10, 1, 0.5, 0.4, 0.1, 10, 1, 0.5)
  totaltime <- 1
  M <- 10
  mainland_n <- M
  island_type <- "oceanic"
  verbose <- FALSE
  replicates <- 1
  sample_freq <- 25
  set.seed(1)
  island_replicates <- list()
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
    formated_CS_sim <- DAISIE:::DAISIE_format_CS(
      island_replicates = island_replicates,
      time = totaltime,
      M = mainland_n,
      sample_freq = sample_freq,
      island_type = island_type,
      verbose = verbose
    )
  )

  expect_equal(
    names(formated_CS_sim[[1]][[1]]),
    c(
      "island_age",
      "not_present_type1",
      "not_present_type2",
      "stt_all",
      "stt_type1",
      "stt_type2"
    )
  )

  expect_equal(
    names(formated_CS_sim[[1]][[2]]),
    c(
      "branching_times",
      "stac",
      "missing_species",
      "other_clades_same_ancestor",
      "type1or2"
    )
  )

  expect_equal(
    length(formated_CS_sim[[1]]), 8
  )

  # Sampled STT has the correct size
  expect_equal(
    dim(formated_CS_sim[[1]][[1]]$stt_all), c(26, 5)
  )
  expect_equal(
    dim(formated_CS_sim[[1]][[1]]$stt_type1), c(26, 5)
  )
  expect_equal(
    dim(formated_CS_sim[[1]][[1]]$stt_type2), c(26, 5)
  )
})

test_that("sampled stt, 1 type, no geodynamics, nonoceanic (same arguments as geodynamics, 5 pars)", {
  totaltime <- 5
  island_age <- 0.4
  clado_rate <- 2.550687345 # cladogenesis rate
  ext_rate <- 2.683454548 # extinction rate
  clade_carr_cap <- 10.0  # clade-level carrying capacity
  imm_rate <- 0.00933207 # immigration rate
  ana_rate <- 1.010073119 # anagenesis rate
  ddmodel_sim <- 11
  pars <- c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate)
  island_type <- "nonoceanic"
  nonoceanic_pars <- c(0.1, 0.9)
  sample_freq <- 25
  mainland_n <- 1
  verbose <- FALSE

  set.seed(1)
  island_replicates <- list()
  out <- list()
  out[[1]] <- DAISIE:::DAISIE_sim_core(
    time = totaltime,
    pars = pars,
    ddmodel_sim = ddmodel_sim,
    mainland_n = mainland_n,
    nonoceanic_pars = nonoceanic_pars,
    island_type = island_type,
  )
  island_replicates[[1]] <- out
  expect_silent(
    formated_CS_sim <- DAISIE:::DAISIE_format_CS(
      island_replicates = island_replicates,
      time = totaltime,
      M = mainland_n,
      sample_freq = sample_freq,
      island_type = island_type,
      verbose = verbose
    )
  )

})

test_that("sampled stt, 2 type, no geodynamics, nonoceanic (same arguments as geodynamics, 10 pars)", {
  skip("DAISIE_sim_min_type2 can't run with nonoceanic")
})

test_that("sampled stt, 1 type, no geodynamics, oceanic (same arguments as geodynamics, 5 pars) - 2 replicates", {
pars <- c(0.4, 0.2, 10, 2, 0.5)
totaltime <- 5
mainland_n <- 1
verbose <- FALSE
sample_freq <- Inf
island_type <- "oceanic"
set.seed(1)
island_replicates <- list()
out <- list()
out[[1]] <- DAISIE:::DAISIE_sim_core(
  time = totaltime,
  pars = pars,
  mainland_n = mainland_n
)
out[[2]] <- DAISIE:::DAISIE_sim_core(
  time = totaltime,
  pars = pars,
  mainland_n = mainland_n
)
island_replicates <- out
expect_silent(
  formated_CS_sim <- DAISIE:::DAISIE_format_CS_full_stt(
    island_replicates = island_replicates,
    time = totaltime,
    M = mainland_n,
    sample_freq = sample_freq,
    island_type = island_type,
    verbose = verbose
  )
)
})



test_that("complete stt, 1 type, no geodynamics, oceanic island (same arguments as geodynamics, 5 pars)", {
  pars <- c(0.4, 0.2, 10, 2, 0.5)
  totaltime <- 1
  mainland_n <- 1
  verbose <- FALSE
  sample_freq <- Inf
  island_type <- "oceanic"
  set.seed(1)

  island_replicates <- list()
  out <- list()
  out[[1]] <- DAISIE:::DAISIE_sim_core(
    time = totaltime,
    pars = pars,
    mainland_n = mainland_n
  )
  island_replicates[[1]] <- out
  expect_silent(
    formated_CS_sim <- DAISIE:::DAISIE_format_CS_full_stt(
      island_replicates = island_replicates,
      time = totaltime,
      M = mainland_n,
      sample_freq = sample_freq,
      island_type = island_type,
      verbose = verbose
    )
  )
  expect_equal(
    formated_CS_sim[[1]][[1]]$island_age,
    5
  )
  expect_equal(
    formated_CS_sim[[1]][[1]]$not_present,
    0
  )
  expect_equal(
    formated_CS_sim[[1]][[1]]$stt_all[2, ],
    c(Time = 0.62240908343582735, nI = 1.0, nA = 0.0, nC = 0.0, present = 1.0)
  )
  expect_equal(
    formated_CS_sim[[1]][[1]]$stt_all[5, ],
    c(Time = 3.8154825768724887, nI = 0.0, nA = 1.0, nC = 0.0, present = 1.0)
  )
  expect_equal(
    formated_CS_sim[[1]][[1]]$stt_all[19, ],
    c(Time = 0.09210138119067679, nI = 1.0, nA = 1.0, nC = 2.0, present = 1.0)
  )

  expect_equal(
    formated_CS_sim[[1]][[2]]$branching_times,
    c(5.0, 1.348741816972570007, 0.092101381190680301)
  )

  expect_equal(
    formated_CS_sim[[1]][[2]]$stac,
    3
  )

  expect_equal(
    formated_CS_sim[[1]][[2]]$missing_species,
    0
  )
})
test_that("complete stt, 1 type, geodynamics, oceanic island (same arguments as no geodynamics, 5 pars)", {

})
test_that("complete stt, 2 type, no geodynamics, oceanic island (same arguments as geodynamics, 10 pars)", {


})
test_that("complete stt, 2 type, geodynamics, oceanic island(same arguments as geodynamics, 10 pars)", {

})
test_that("complete stt, 1 type, no geodynamics, nonoceanic (same arguments as geodynamics, 5 pars)", {

})
test_that("complete stt, 2 type, no geodynamics, nonoceanic (same arguments as geodynamics, 10 pars)", {

})


