context("DAISIE_format_IW")

test_that("silent with empty island with correct output", {
  pars <- c(0.4, 0.2, 10, 0.0001, 0.5)
  time <- 1
  mainland_n <- 10
  verbose <- FALSE
  sample_freq <- 1
  start_midway <- FALSE
  set.seed(1)
  island_replicates <- list()
  island_replicates[[1]] <- DAISIE:::DAISIE_sim_core_constant_rate(
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
  expected_IW_format <- list()
  expected_IW_format[[1]] <- list()
  stt_all <- matrix(ncol = 4, nrow = 2)
  colnames(stt_all) <- c("Time", "nI", "nA", "nC")
  stt_all[1, ] <- c(1, 0, 0, 0)
  stt_all[2, ] <- c(0, 0, 0, 0)
  brts_table <- matrix(ncol = 4, nrow = 1)
  colnames(brts_table) <- c("brt", "clade", "event", "endemic")
  brts_table[1, ] <- c(1, 0, 0, NA)
  expected_IW_format[[1]][[1]] <- list(island_age = 1,
                                       not_present = 10,
                                       stt_all = stt_all,
                                       brts_table = brts_table)
  expect_true(all.equal(formated_IW_sim, expected_IW_format, tolerance = 1e-7))
})

test_that("silent with non-empty island with correct output", {
  pars <- c(0.4, 0.2, 10, 1, 0.5)
  time <- 1
  mainland_n <- 10
  verbose <- FALSE
  sample_freq <- 1
  start_midway <- FALSE
  set.seed(1)
  island_replicates <- list()
  island_replicates[[1]] <- DAISIE:::DAISIE_sim_core_constant_rate(
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
  expected_IW_format <- list()
  expected_IW_format[[1]] <- list()
  stt_all <- matrix(ncol = 4, nrow = 2)
  colnames(stt_all) <- c("Time", "nI", "nA", "nC")
  stt_all[1, ] <- c(1, 0, 0, 0)
  stt_all[2, ] <- c(0, 2, 0, 3)
  brts_table <- matrix(ncol = 4, nrow = 6)
  colnames(brts_table) <- c("brt", "clade", "event", "endemic")
  brts_table[1, ] <- c(1, 0, 0, NA)
  brts_table[2, ] <- c(0.9244818166871660, 1, 1, 1)
  brts_table[3, ] <- c(0.9105856673960619, 1, 2, 1)
  brts_table[4, ] <- c(0.5557734125062590, 2, 1, 0)
  brts_table[5, ] <- c(0.5288428248966160, 3, 1, 0)
  brts_table[6, ] <- c(0.3146835586399670, 1, 3, 1)
  expected_IW_format[[1]][[1]] <- list(island_age = 1,
                                       not_present = 3,
                                       stt_all = stt_all,
                                       brts_table = brts_table)
  expected_IW_format[[1]][[2]] <- list(branching_times = c(1.00000000000000,
                                                           0.924481816687166,
                                                           0.910585667396062,
                                                           0.314683558639967),
                                       stac = 2,
                                       missing_species = 0)
  expected_IW_format[[1]][[3]] <- list(branching_times = c(1.000000000000000,
                                                           0.555773412506259),
                                       stac = 4,
                                       missing_species = 0)
  expected_IW_format[[1]][[4]] <- list(branching_times = c(1.000000000000000,
                                                           0.5288428248966160),
                                       stac = 4,
                                       missing_species = 0)
  expect_true(all.equal(formated_IW_sim, expected_IW_format, tolerance = 1e-7))
})

test_that("DAISIE_format_IW prints when verbose = TRUE", {
  pars <- c(0.4, 0.2, 10, 1, 0.5)
  time <- 1
  mainland_n <- 1000
  verbose <- TRUE
  sample_freq <- 1
  start_midway <- FALSE
  set.seed(1)
  island_replicates <- list()
  island_replicates[[1]] <- DAISIE:::DAISIE_sim_core_constant_rate(
    time = time,
    pars = pars,
    mainland_n = mainland_n
  )
  expect_output(
    formated_IW_sim <- DAISIE:::DAISIE_format_IW(
      island_replicates = island_replicates,
      time = time,
      M = mainland_n,
      sample_freq = sample_freq,
      verbose = verbose
    ),
    "Island being formatted: 1/1"
  )
})

test_that("silent with empty nonoceanic island with correct output", {
  pars <- c(0.4, 2, 10, 0.0001, 0.5)
  time <- 1
  mainland_n <- 10
  nonoceanic_pars <- c(0.2, 0.5)
  verbose <- FALSE
  sample_freq <- 1
  start_midway <- FALSE
  set.seed(1)
  island_replicates <- list()
  island_replicates[[1]] <- DAISIE:::DAISIE_sim_core_constant_rate(
    time = time,
    mainland_n = mainland_n,
    pars = pars,
    nonoceanic_pars = nonoceanic_pars
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
  expected_IW_format <- list()
  expected_IW_format[[1]] <- list()
  stt_all <- matrix(ncol = 4, nrow = 2)
  colnames(stt_all) <- c("Time", "nI", "nA", "nC")
  stt_all[1, ] <- c(1, 1, 2, 0)
  stt_all[2, ] <- c(0, 0, 0, 0)
  brts_table <- matrix(ncol = 4, nrow = 1)
  colnames(brts_table) <- c("brt", "clade", "event", "endemic")
  brts_table[1, ] <- c(1, 0, 0, NA)
  expected_IW_format[[1]][[1]] <- list(island_age = 1,
                                       not_present = 10,
                                       stt_all = stt_all,
                                       brts_table = brts_table)
  expect_equal(formated_IW_sim, expected_IW_format)
})

test_that("silent with non-empty nonoceanic island with
          correct output", {
            pars <- c(0.4, 0.2, 1, 1, 0.5)
            time <- 1
            mainland_n <- 10
            nonoceanic_pars <- c(0.2, 0.5)
            verbose <- FALSE
            sample_freq <- 1
            start_midway <- FALSE
            set.seed(1)
            island_replicates <- list()
            island_replicates[[1]] <- DAISIE:::DAISIE_sim_core_constant_rate(
              time = time,
              mainland_n = mainland_n,
              pars = pars,
              nonoceanic_pars = nonoceanic_pars
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
            expected_IW_format <- list()
            expected_IW_format[[1]] <- list()
            stt_all <- matrix(ncol = 4, nrow = 2)
            colnames(stt_all) <- c("Time", "nI", "nA", "nC")
            stt_all[1, ] <- c(1, 1, 2, 0)
            stt_all[2, ] <- c(0, 0, 2, 0)
            brts_table <- matrix(ncol = 4, nrow = 3)
            colnames(brts_table) <- c("brt", "clade", "event", "endemic")
            brts_table[1, ] <- c(1, 0, 0, NA)
            brts_table[2, ] <- c(1, 2, 1, 1)
            brts_table[3, ] <- c(1, 1, 1, 1)
            expected_IW_format[[1]][[1]] <- list(island_age = 1,
                                                 not_present = 2,
                                                 stt_all = stt_all,
                                                 brts_table = brts_table)
            expected_IW_format[[1]][[2]] <- list(branching_times = c(1,
                                                                     1),
                                                 stac = 2,
                                       missing_species = 0)
  expected_IW_format[[1]][[3]] <- list(branching_times = c(1,
                                                           1),
                                       stac = 2,
                                       missing_species = 0)
  expect_true(all.equal(formated_IW_sim, expected_IW_format, tolerance = 1e-7))
})

test_that("Add_brt_table output is correct when length(island) == 1", {
  stt_all <- matrix(ncol = 4, nrow = 2)
  colnames(stt_all) <- c("Time", "nI", "nA", "nC")
  stt_all[1, ] <- c(1, 0, 0, 0)
  stt_all[2, ] <- c(0, 0, 0, 0)
  island <- list()
  island[[1]] <- list(island_age = 1,
                 not_present = 100,
                 stt_all = stt_all,
                 init_nonend_spec = 0,
                 init_end_spec = 0)
  formatted_brt <- DAISIE:::Add_brt_table(island)
  brt_table <- matrix(ncol = 4, nrow = 1)
  colnames(brt_table) <- c("brt", "clade", "event", "endemic")
  brt_table[1, ] <- c(1, 0, 0, NA)
  expected_brt <- list()
  expected_brt[[1]] <- list(island_age = 1,
                        not_present = 100,
                        stt_all = stt_all,
                        init_nonend_spec = 0,
                        init_end_spec = 0,
                        brts_table = brt_table)
  expect_true(all.equal(formatted_brt, expected_brt))
})
test_that("Add_brt_table output is correct when length(island) != 1", {
  stt_all <- matrix(ncol = 4, nrow = 2)
  colnames(stt_all) <- c("Time", "nI", "nA", "nC")
  stt_all[1, ] <- c(1, 0, 0, 0)
  stt_all[2, ] <- c(0, 2, 0, 3)
  island <- list()
  island[[1]] <- list(island_age = 1,
                      not_present = 3,
                      stt_all = stt_all,
                      init_nonend_spec = 0,
                      init_end_spec = 0)
  island[[2]] <- list(branching_times = c(1.0000000,
                                          0.5557734),
                      stac = 4,
                      missing_species = 0)
  island[[3]] <- list(branching_times = c(1.0000000,
                                          0.9244818,
                                          0.9105857,
                                          0.3146836),
                      stac = 2,
                      missing_species = 0)
  island[[4]] <- list(brancing_times = c(1.0000000,
                                         0.5288428),
                      stac = 4,
                      missing_species = 0)
  formatted_brt <- DAISIE:::Add_brt_table(island)
  brt_table <- matrix(ncol = 4, nrow = 5)
  colnames(brt_table) <- c("brt", "clade", "event", "endemic")
  brt_table[1, ] <- c(1, 0, 0, NA)
  brt_table[2, ] <- c(0.9244818, 1, 1, 1)
  brt_table[3, ] <- c(0.9105857, 1, 2, 1)
  brt_table[4, ] <- c(0.5557734, 2, 1, 0)
  brt_table[5, ] <- c(0.3146836, 1, 3, 1)
  expected_brt <- list()
  expected_brt[[1]] <- list(island_age = 1,
                        not_present = 3,
                        stt_all = stt_all,
                        init_nonend_spec = 0,
                        init_end_spec = 0,
                        brts_table = brt_table)
  expected_brt[[2]] <- list(branching_times = c(1.0000000,
                                                 0.9244818,
                                                 0.9105857,
                                                 0.3146836),
                             stac = 2,
                             missing_species = 0)
  expected_brt[[3]] <- list(branching_times = c(1.0000000,
                                                 0.5557734),
                             stac = 4,
                             missing_species = 0)
  expect_equal(formatted_brt, expected_brt)
})
#test_that("Add_brt_table output is correct when length(stac1_5s) != 0")
#test_that("Add_brt_table output is correct when length(stac1_5s) == 0")
#test_that("Add_brt_table output is correct when length(island_no_stac1or5) != 0")

test_that("abuse", {
  expect_error(DAISIE:::DAISIE_format_IW("nonsense"))
})

test_that("abuse", {
  expect_error(DAISIE:::Add_brt_table("nonsense"))
})


######  add format_IW_trait tests from here ####
test_that("silent with empty island with correct output", {
  pars <- c(0.4, 0.2, 10, 0.0001, 0.5)
  trait_pars <- create_trait_pars(
    trans_rate = 0,
    immig_rate2 = 0.0002,
    ext_rate2 = 0.2,
    ana_rate2 = 0.5,
    clado_rate2 = 0.4,
    trans_rate2 = 0,
    M2 = 10)
  time <- 1
  mainland_n <- 10
  verbose <- FALSE
  sample_freq <- 1
  start_midway <- FALSE
  set.seed(1)
  island_replicates <- list()
  island_replicates[[1]] <- DAISIE:::DAISIE_sim_core_trait_dependent(
    time = time,
    pars = pars,
    trait_pars = trait_pars,
    mainland_n = mainland_n
  )
  expect_silent(
    formated_IW_sim <- DAISIE:::DAISIE_format_IW(
      island_replicates = island_replicates,
      time = time,
      M = mainland_n,
      sample_freq = sample_freq,
      verbose = verbose,
      trait_pars = trait_pars
    )
  )
  expected_IW_format <- list()
  expected_IW_format[[1]] <- list()
  stt_all <- matrix(ncol = 7, nrow = 2)
  colnames(stt_all) <- c("Time", "nI", "nA", "nC", "nI2", "nA2", "nC2")
  stt_all[1, ] <- c(1, 0, 0, 0, 0, 0, 0)
  stt_all[2, ] <- c(0, 0, 0, 0, 0, 0, 0)
  brts_table <- matrix(ncol = 4, nrow = 1)
  colnames(brts_table) <- c("brt", "clade", "event", "endemic")
  brts_table[1, ] <- c(1, 0, 0, NA)
  expected_IW_format[[1]][[1]] <- list(island_age = 1,
                                       not_present = 20,
                                       stt_all = stt_all,
                                       brts_table = brts_table)
  expect_true(all.equal(formated_IW_sim, expected_IW_format, tolerance = 1e-7))
})

test_that("silent when species with two trait states with
          correct output", {
            pars <- c(0.4, 0.2, 10, 0.06, 0.5)
            time <- 5
            mainland_n <- 10
            nonoceanic_pars <- c(0, 0)
            verbose <- FALSE
            replicates <- 3
            island_ontogeny = 0
            sea_level = 0
            extcutoff = 1000
            sample_freq <- 1
            trait_pars <- create_trait_pars(
              trans_rate = 0,
              immig_rate2 = 0.1,
              ext_rate2 = 0.2,
              ana_rate2 = 0.5,
              clado_rate2 = 0.4,
              trans_rate2 = 0,
              M2 = 5)
            island_replicates <- list()
            verbose <- FALSE
            set.seed(1)
            island_replicates[[1]] <- DAISIE:::DAISIE_sim_core_trait_dependent(
              time = time,
              mainland_n = mainland_n,
              pars = pars,
              nonoceanic_pars = nonoceanic_pars,
              trait_pars = trait_pars,
              island_ontogeny = island_ontogeny,
              sea_level = sea_level,
              extcutoff = extcutoff
            )
            expect_silent(
              formated_IW_sim <- DAISIE:::DAISIE_format_IW(
                island_replicates = island_replicates,
                time = time,
                M = mainland_n,
                sample_freq = sample_freq,
                verbose = verbose,
                trait_pars = trait_pars
              )
            )
            expected_IW_format <- list()
            expected_IW_format[[1]] <- list()
            stt_all <- matrix(ncol = 7, nrow = 2)
            colnames(stt_all) <- c("Time", "nI", "nA", "nC", "nI2", "nA2", "nC2")
            stt_all[1, ] <- c(5, 0, 0, 0, 0, 0, 0)
            stt_all[2, ] <- c(0, 0, 1, 2, 0, 0, 0)
            brts_table <- matrix(ncol = 4, nrow = 4)
            colnames(brts_table) <- c("brt", "clade", "event", "endemic")
            brts_table[1, ] <- c(5.000000000, 0, 0, NA)
            brts_table[2, ] <- c(3.102613675, 1, 1, 1)
            brts_table[3, ] <- c(1.505629998, 2, 1, 1)
            brts_table[4, ] <- c(1.262456559, 2, 2, 1)
            expected_IW_format[[1]][[1]] <- list(island_age = 5,
                                                 not_present = 13,
                                                 stt_all = stt_all,
                                                 brts_table = brts_table)
            expected_IW_format[[1]][[2]] <- list(branching_times = c(5.000000000,
                                                                     3.102613675),
                                                 stac = 2,
                                                 missing_species = 0)
            expected_IW_format[[1]][[3]] <- list(branching_times = c(5.000000000,
                                                                     1.505629998,
                                                                     1.262456559),
                                                 stac = 2,
                                                 missing_species = 0)
            expect_true(all.equal(formated_IW_sim, expected_IW_format, tolerance = 1e-7))
          })
