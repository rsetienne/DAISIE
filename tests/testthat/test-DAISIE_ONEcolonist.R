context("DAISIE_ONEcolonist")

test_that("DAISIE_ONEcolonist works on an oceanic DAISIE_sim_core", {
  sim_time <- 10
  n_mainland_species <- 1
  clado_rate <- 1.0
  ext_rate <- 0.1
  carr_cap <- 4
  imm_rate <- 1.0
  ana_rate <- 1.0
  set.seed(1)
  area_pars <- DAISIE::create_area_pars(
    max_area = 1,
    current_area = 1,
    proportional_peak_t = 0,
    total_island_age = 0,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0)
  hyper_pars <- create_hyper_pars(d = 0, x = 0)
  nonoceanic_pars <- c(0, 0)
  sim <- DAISIE:::DAISIE_sim_core_constant_rate(
    time = sim_time,
    mainland_n = n_mainland_species,
    pars = c(clado_rate, ext_rate, carr_cap, imm_rate, ana_rate),
    area_pars = area_pars,
    hyper_pars = hyper_pars,
    nonoceanic_pars = nonoceanic_pars)
  stt_table <- sim$stt_table
  island_spec <- matrix(nrow = 4, ncol = 7, data = "x")
  island_spec[, 1] <- c("6", "10", "9", "11")
  island_spec[, 2] <- c("1", "1", "1", "1")
  island_spec[, 3] <- c("0.755181833128345",
                        "0.755181833128345",
                        "0.755181833128345",
                        "0.755181833128345")
  island_spec[, 4] <- c("C", "C", "C", "C")
  island_spec[, 5] <- c("AA", "ABA", "B", "ABB")
  island_spec[, 6] <- c("0.755181833128345",
                        "2.66196121187029",
                        "0.808949241539306",
                        "9.7444563652974")
  island_spec[, 7] <- c(NA, NA, NA, NA)
  colnames(island_spec) <- c(
    "Species",
    "Mainland Ancestor",
    "Colonisation time (BP)",
    "Species type",
    "branch_code",
    "branching time (BP)",
    "Anagenetic_origin"
  )
  init_nonend_spec <- sim$init_nonend_spec
  init_end_spec <- sim$init_end_spec
  carrying_capacity <- sim$carrying_capacity
  result <- DAISIE:::DAISIE_ONEcolonist(
    time = sim_time,
    island_spec = island_spec,
    stt_table = stt_table
  )
  expect_equal(result$stt_table, stt_table)
  expect_true(
    all.equal(
      result$branching_times,
      c(10.0000000, 9.7444564, 2.6619612, 0.8089492, 0.7551818),
      tolerance = 1.0e-7
    )
  )
  expect_equal(result$stac, sim$stac)
  expect_equal(result$missing_species, sim$missing_species)
})

#test_that("DAISIE_ONEcolonist works on a nonoceanic DAISIE_sim_core")
#test_that("DAISIE_ONEcolonist works on an oceanic DAISIE_sim_core with
#other_clades_same_ancestor)
#test_that("DAISIE_ONEcolonist works on a nonoceanic DAISIE_sim_core with
#other_clades_same_ancestor)

test_that("DAISIE_ONEcolonist works with >=2 cladogenetic with same ancestor", {
  set.seed(42)
  sim_time <- 10
  n_mainland_species <- 1
  clado_rate <- 1
  ext_rate <- 0.00001
  carr_cap <- 4
  imm_rate <- 1
  ana_rate <- 0.000001
  area_pars <- DAISIE::create_area_pars(
    max_area = 1,
    current_area = 1,
    proportional_peak_t = 0,
    total_island_age = 0,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0)
  hyper_pars <- create_hyper_pars(d = 0, x = 0)
  nonoceanic_pars <- c(0, 0)
  expect_silent(out <- DAISIE:::DAISIE_sim_core_constant_rate(
    time = sim_time,
    mainland_n = n_mainland_species,
    pars = c(clado_rate, ext_rate, carr_cap, imm_rate, ana_rate),
    area_pars = area_pars,
    hyper_pars = hyper_pars,
    nonoceanic_pars = nonoceanic_pars
  )
  )
})


test_that("DAISIE_ONEcolonist works with >=2 anagenetic with same ancestor", {
  set.seed(42)
  sim_time <- 10
  n_mainland_species <- 1
  clado_rate <- 0.0000001
  ext_rate <- 0.00001
  carr_cap <- 4
  imm_rate <- 1
  ana_rate <- 2
  area_pars <- DAISIE::create_area_pars(
    max_area = 1,
    current_area = 1,
    proportional_peak_t = 0,
    total_island_age = 0,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0)
  hyper_pars <- create_hyper_pars(d = 0, x = 0)
  nonoceanic_pars <- c(0, 0)
  expect_silent(out <- DAISIE:::DAISIE_sim_core_constant_rate(
    time = sim_time,
    mainland_n = n_mainland_species,
    pars = c(clado_rate, ext_rate, carr_cap, imm_rate, ana_rate),
    area_pars = area_pars,
    hyper_pars = hyper_pars,
    nonoceanic_pars = nonoceanic_pars)
  )
})

test_that("DAISIE_ONEcolonist works with >=2 nonendemic with same ancestor", {
  set.seed(44)
  sim_time <- 10
  n_mainland_species <- 1
  clado_rate <- 0.0000001
  ext_rate <- 0.00001
  carr_cap <- 4
  imm_rate <- 3
  ana_rate <- 1
  area_pars <- DAISIE::create_area_pars(
    max_area = 1,
    current_area = 1,
    proportional_peak_t = 0,
    total_island_age = 0,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0)
  hyper_pars <- create_hyper_pars(d = 0, x = 0)
  nonoceanic_pars <- c(0, 0)
  expect_silent(out <- DAISIE:::DAISIE_sim_core_constant_rate(
    time = sim_time,
    mainland_n = n_mainland_species,
    pars = c(clado_rate, ext_rate, carr_cap, imm_rate, ana_rate),
    area_pars = area_pars,
    hyper_pars = hyper_pars,
    nonoceanic_pars = nonoceanic_pars
  )
  )
})

test_that("DAISIE_ONEcolonist stac and brts works for single colonist", {
  sim_time <- 10

  island_spec <- matrix(nrow = 4, ncol = 7, data = "x")
  island_spec[, 1] <- c("6", "10", "9", "11")
  island_spec[, 2] <- c("1", "1", "1", "1")
  island_spec[, 3] <- c("0.755181833128345",
                        "0.755181833128345",
                        "0.755181833128345",
                        "0.755181833128345")
  island_spec[, 4] <- c("C", "C", "C", "C")
  island_spec[, 5] <- c("AA", "ABA", "B", "ABB")
  island_spec[, 6] <- c("0.755181833128345",
                        "2.66196121187029",
                        "0.808949241539306",
                        "9.7444563652974")
  island_spec[, 7] <- c(NA, NA, NA, NA)
  colnames(island_spec) <- c(
    "Species",
    "Mainland Ancestor",
    "Colonisation time (BP)",
    "Species type",
    "branch_code",
    "branching time (BP)",
    "Anagenetic_origin"
  )
  stt_table <- NULL
  result <- DAISIE:::DAISIE_ONEcolonist(
    time = sim_time,
    island_spec = island_spec,
    stt_table = stt_table
  )
  expect_equal(
    result$branching_times,
    c(sim_time, sort(as.numeric(island_spec[, 6]), decreasing = TRUE))
  )
  # Stac 2 for regular clade on island with no recolonization
  expect_equal(result$stac, 2)
})

test_that("DAISIE_ONEcolonist stac and brts works for 1 nonendemic colonist", {
  # With just 1 nonendemic recolonist, function works

  sim_time <- 2

  island_spec <- matrix(nrow = 2, ncol = 7, data = "x")
  island_spec[, 1] <- c("1", "2")
  island_spec[, 2] <- c("1", "1")
  island_spec[, 3] <- c("0.7", "0.6")
  island_spec[, 4] <- c("A", "I")
  island_spec[, 5] <- c(NA, NA)
  island_spec[, 6] <- c(NA, NA)
  island_spec[, 7] <- c("Immig_parent", NA)
  colnames(island_spec) <- c(
    "Species",
    "Mainland Ancestor",
    "Colonisation time (BP)",
    "Species type",
    "branch_code",
    "branching time (BP)",
    "Anagenetic_origin"
  )
  stt_table <- NULL
  result <- DAISIE:::DAISIE_ONEcolonist(
    time = sim_time,
    island_spec = island_spec,
    stt_table = stt_table
  )
  expect_equal(
    result$branching_times,
    c(sim_time, max(as.numeric(island_spec[, 3])))
  )

  # Stac 3 for recolonisation cases
  expect_equal(result$stac, 3)

  expect_equal(result$all_colonisations[[1]]$event_times, c(2, 0.7))
  expect_equal(result$all_colonisations[[1]]$species_type, "A")

  expect_equal(result$all_colonisations[[2]]$event_times, c(2, 0.6))
  expect_equal(result$all_colonisations[[2]]$species_type, "I")

})

test_that("DAISIE_ONEcolonist stac and brts works for 2 endemic colonists", {
  # With > 1 endemic recolonist, function works

  sim_time <- 2

  island_spec <- matrix(nrow = 3, ncol = 7, data = "x")
  island_spec[, 1] <- c("1", "2", "3")
  island_spec[, 2] <- c("1", "1", "1")
  island_spec[, 3] <- c("0.7", "0.6", "0.5")
  island_spec[, 4] <- c("A", "A", "I")
  island_spec[, 5] <- c(NA, NA, NA)
  island_spec[, 6] <- c(NA, NA, NA)
  island_spec[, 7] <- c("Immig_parent", "Immig_parent", NA)
  colnames(island_spec) <- c(
    "Species",
    "Mainland Ancestor",
    "Colonisation time (BP)",
    "Species type",
    "branch_code",
    "branching time (BP)",
    "Anagenetic_origin"
  )
  stt_table <- NULL
  result <- DAISIE:::DAISIE_ONEcolonist(
    time = sim_time,
    island_spec = island_spec,
    stt_table = stt_table
  )
  # Only include oldest colonisation time time
  expect_equal(
    result$branching_times,
    c(sim_time, as.numeric(island_spec[, 3])[1:2])
  )
  # stac 3 for recolonisation cases
  expect_equal(result$stac, 3)

  expect_equal(result$all_colonisations[[1]]$event_times, c(2, 0.7))
  expect_equal(result$all_colonisations[[1]]$species_type, "A")

  expect_equal(result$all_colonisations[[2]]$event_times, c(2, 0.6))
  expect_equal(result$all_colonisations[[2]]$species_type, "A")

  expect_equal(result$all_colonisations[[3]]$event_times, c(2, 0.5))
  expect_equal(result$all_colonisations[[3]]$species_type, "I")
})

test_that("DAISIE_ONEcolonist stac and brts works for 2 endemic clades", {
  # With > 1 endemic clades, function works

  sim_time <- 2


  # Species Mainland Ancestor Colonisation time (BP) Species type branch_code branching time (BP) Anagenetic_origin
  # [1,] "4"     "1"               "1.13468671408026"     "C"          "AA"        "1.13468671408026"  NA
  # [2,] "3"     "1"               "1.13468671408026"     "C"          "B"         "0.96545899791969"  NA
  # [3,] "5"     "1"               "1.13468671408026"     "C"          "AB"        "0.68696590746724"  NA
  # [4,] "6"     "1"               "0.67395467208331"     "C"          "A"         "0.67395467208331"  NA
  # [5,] "7"     "1"               "0.67395467208331"     "C"          "B"         "0.34198900695798"  NA


  island_spec <- matrix(nrow = 5, ncol = 7, data = "x")
  island_spec[, 1] <- c("4", "3", "5", "6", "7")
  island_spec[, 2] <- c("1", "1", "1", "1", "1")
  island_spec[, 3] <-
    c("1.13468671408026",
      "1.13468671408026",
      "1.13468671408026",
      "0.67395467208331",
      "0.67395467208331")
  island_spec[, 4] <- c("C", "C", "C", "C", "C")
  island_spec[, 5] <- c("AA", "B", "AB", "A", "B")
  island_spec[, 6] <-
    c(1.13468671408026,
      0.96545899791969,
      0.68696590746724,
      0.67395467208331,
      0.34198900695798)
  island_spec[, 7] <- c(NA, NA, NA, NA, NA)
  colnames(island_spec) <- c(
    "Species",
    "Mainland Ancestor",
    "Colonisation time (BP)",
    "Species type",
    "branch_code",
    "branching time (BP)",
    "Anagenetic_origin"
  )
  stt_table <- NULL
  result <- DAISIE:::DAISIE_ONEcolonist(
    time = sim_time,
    island_spec = island_spec,
    stt_table = stt_table
  )
  btimes <- sort(as.numeric(island_spec[, 6]), decreasing = TRUE)
  btimes_sans_yng_col <- btimes[-4]
  expect_equal(
    result$branching_times,
    c(sim_time, btimes_sans_yng_col)
  )
  # stac 3 for recolonisation cases
  expect_equal(result$stac, 3)

  # all_colonisations
  expect_equal(result$all_colonisations[[1]]$event_times, c(
    2.0,
    1.13468671408026,
    0.96545899791969,
    0.68696590746724
  ))
  expect_equal(result$all_colonisations[[1]]$species_type, "C")

  expect_equal(result$all_colonisations[[2]]$event_times, c(
    2.0,
    0.67395467208331,
    0.34198900695798
  ))
  expect_equal(result$all_colonisations[[2]]$species_type, "C")
})

test_that("DAISIE_ONEcolonist stac and brts works for 2 endemic clades", {
  # With > 1 endemic clades, function works

  sim_time <- 2


  # Species Mainland Ancestor Colonisation time (BP) Species type branch_code branching time (BP) Anagenetic_origin
  # [1,] "4"     "1"               "1.13468671408026"     "C"          "AA"        "1.13468671408026"  NA
  # [2,] "3"     "1"               "1.13468671408026"     "C"          "B"         "0.96545899791969"  NA
  # [3,] "5"     "1"               "1.13468671408026"     "C"          "AB"        "0.68696590746724"  NA
  # [4,] "6"     "1"               "0.67395467208331"     "C"          "A"         "0.67395467208331"  NA
  # [5,] "7"     "1"               "0.67395467208331"     "C"          "B"         "0.34198900695798"  NA


  island_spec <- matrix(nrow = 5, ncol = 7, data = "x")
  island_spec[, 1] <- c("4", "3", "5", "6", "7")
  island_spec[, 2] <- c("1", "1", "1", "1", "1")
  island_spec[, 3] <-
    c("1.13468671408026",
      "1.13468671408026",
      "1.13468671408026",
      "0.67395467208331",
      "0.67395467208331")
  island_spec[, 4] <- c("C", "C", "C", "C", "C")
  island_spec[, 5] <- c("AA", "B", "AB", "A", "B")
  island_spec[, 6] <-
    c(1.13468671408026,
      0.96545899791969,
      0.68696590746724,
      0.67395467208331,
      0.34198900695798)
  island_spec[, 7] <- c(NA, NA, NA, NA, NA)
  colnames(island_spec) <- c(
    "Species",
    "Mainland Ancestor",
    "Colonisation time (BP)",
    "Species type",
    "branch_code",
    "branching time (BP)",
    "Anagenetic_origin"
  )
  stt_table <- NULL
  result <- DAISIE:::DAISIE_ONEcolonist(
    time = sim_time,
    island_spec = island_spec,
    stt_table = stt_table
  )
  btimes <- sort(as.numeric(island_spec[, 6]), decreasing = TRUE)
  btimes_sans_yng_col <- btimes[-4]
  expect_equal(
    result$branching_times,
    c(sim_time, btimes_sans_yng_col)
  )
  # stac 3 for recolonisation cases
  expect_equal(result$stac, 3)

  # all_colonisations
  expect_equal(result$all_colonisations[[1]]$event_times, c(
    2.0,
    1.13468671408026,
    0.96545899791969,
    0.68696590746724
  ))
  expect_equal(result$all_colonisations[[1]]$species_type, "C")

  expect_equal(result$all_colonisations[[2]]$event_times, c(
    2.0,
    0.67395467208331,
    0.34198900695798
  ))
  expect_equal(result$all_colonisations[[2]]$species_type, "C")
})

