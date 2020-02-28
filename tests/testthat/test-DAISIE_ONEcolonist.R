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
  sim <- DAISIE:::DAISIE_sim_core_constant_rate(
    time = sim_time,
    mainland_n = n_mainland_species,
    pars = c(clado_rate, ext_rate, carr_cap, imm_rate, ana_rate))
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
  expect_silent(out <- DAISIE:::DAISIE_sim_core_constant_rate(
    time = sim_time,
    mainland_n = n_mainland_species,
    pars = c(clado_rate, ext_rate, carr_cap, imm_rate, ana_rate)
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
  expect_silent(out <- DAISIE:::DAISIE_sim_core_constant_rate(
    time = sim_time,
    mainland_n = n_mainland_species,
    pars = c(clado_rate, ext_rate, carr_cap, imm_rate, ana_rate)
  )
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
  expect_silent(out <- DAISIE:::DAISIE_sim_core_constant_rate(
    time = sim_time,
    mainland_n = n_mainland_species,
    pars = c(clado_rate, ext_rate, carr_cap, imm_rate, ana_rate)
  )
  )
})

