context("DAISIE_sim_constant_rate")

test_that("A divdepmodel = 'CS' run should produce no output", {
  n_mainland_species <- 100
  island_age <- 5
  clado_rate <- 1.0
  ext_rate <- 1.0
  clade_carr_cap <- 10.0
  imm_rate <- 0.01
  ana_rate <- 1.0
  expect_silent(
    DAISIE_sim_constant_rate(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
})

test_that("A divdepmodel = 'IW' run should produce no output", {
  n_mainland_species <- 100
  island_age <- 0.4
  clado_rate <- 1.0
  ext_rate <- 1.0
  clade_carr_cap <- 10.0
  imm_rate <- 0.01
  ana_rate <- 1.0
  expect_silent(
    DAISIE_sim_constant_rate(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      divdepmodel = "IW",
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
})

test_that("A divdepmodel = 'GW' run should produce no output", {
  n_mainland_species <- 100
  island_age <- 0.4
  clado_rate <- 1.0
  ext_rate <- 1.0
  clade_carr_cap <- 10.0
  imm_rate <- 0.01
  ana_rate <- 1.0
  num_guilds <- 5
  expect_silent(
    DAISIE_sim_constant_rate(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      divdepmodel = "GW",
      num_guilds = num_guilds,
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
})

test_that("A 2 type with replicates_apply_type2 == TRUE
          divdepmodel = 'CS' run should produce no output", {
            n_mainland_species <- 100
            island_age <- 0.4
            clado_rate_type_1 <- 1.0
            ext_rate_type_1 <- 1.0
            clade_carr_cap_type_1 <- 10.0
            imm_rate_type_1 <- 0.01
            ana_rate_type_1 <- 1.0
            clado_rate_type_2 <- 1.0
            ext_rate_type_2 <- 1.0
            clade_carr_cap_type_2 <- 10.0
            imm_rate_type_2 <- 0.01
            ana_rate_type_2 <- 1.0
            prop_type2_pool <- 0.1
            replicates_apply_type2 <- TRUE
            expect_silent(
              sim <- DAISIE_sim_constant_rate(
                time = island_age,
                M = n_mainland_species,
                pars = c(clado_rate_type_1,
                         ext_rate_type_1,
                         clade_carr_cap_type_1,
                         imm_rate_type_1,
                         ana_rate_type_1,
                         clado_rate_type_2,
                         ext_rate_type_2,
                         clade_carr_cap_type_2,
                         imm_rate_type_2,
                         ana_rate_type_2),
                replicates = 1,
                prop_type2_pool = prop_type2_pool,
                replicates_apply_type2 = replicates_apply_type2,
                plot_sims = FALSE,
                verbose = FALSE
              )
            )
          })


test_that("A 2 type with replicates_apply_type2 == FALSE
          divdepmodel = 'CS' run should produce no output", {
            n_mainland_species <- 100
            island_age <- 0.4
            clado_rate_type_1 <- 1.0
            ext_rate_type_1 <- 1.0
            clade_carr_cap_type_1 <- 10.0
            imm_rate_type_1 <- 0.01
            ana_rate_type_1 <- 1.0
            clado_rate_type_2 <- 1.0
            ext_rate_type_2 <- 1.0
            clade_carr_cap_type_2 <- 10.0
            imm_rate_type_2 <- 0.01
            ana_rate_type_2 <- 1.0
            prop_type2_pool <- 0.1
            replicates_apply_type2 <- FALSE
            expect_silent(
              sim <- DAISIE_sim_constant_rate(
                time = island_age,
                M = n_mainland_species,
                pars = c(clado_rate_type_1,
                         ext_rate_type_1,
                         clade_carr_cap_type_1,
                         imm_rate_type_1,
                         ana_rate_type_1,
                         clado_rate_type_2,
                         ext_rate_type_2,
                         clade_carr_cap_type_2,
                         imm_rate_type_2,
                         ana_rate_type_2),
                replicates = 1,
                prop_type2_pool = prop_type2_pool,
                replicates_apply_type2 = replicates_apply_type2,
                plot_sims = FALSE,
                verbose = FALSE
              )
            )
          })

test_that("output is correct for divdepmodl = 'CS'", {
  n_mainland_species <- 100
  island_age <- 0.4
  clado_rate <- 1.0
  ext_rate <- 1.0
  clade_carr_cap <- 10.0
  imm_rate <- 0.01
  ana_rate <- 1.0
  replicates <- 1
  sim <- DAISIE_sim_constant_rate(
    time = island_age,
    M = n_mainland_species,
    pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
    replicates = replicates,
    plot_sims = FALSE,
    verbose = FALSE
  )
  expect_true(is.list(sim))
  expect_true(length(sim) == replicates)
})

test_that("output is correct for divdepmodel = 'IW'", {
  n_mainland_species <- 100
  island_age <- 0.4
  clado_rate <- 1.0
  ext_rate <- 1.0
  clade_carr_cap <- 10.0
  imm_rate <- 0.01
  ana_rate <- 1.0
  replicates <- 1
  sim <- DAISIE_sim_constant_rate(
    time = island_age,
    M = n_mainland_species,
    pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
    replicates = replicates,
    divdepmodel = "IW",
    plot_sims = FALSE,
    verbose = FALSE
  )
  expect_true(is.list(sim))
  expect_true(length(sim) == replicates)
})

test_that("output is correct for divdepmodl = 'GW'", {
  n_mainland_species <- 100
  island_age <- 0.4
  clado_rate <- 1.0
  ext_rate <- 1.0
  clade_carr_cap <- 10.0
  imm_rate <- 0.01
  ana_rate <- 1.0
  replicates <- 1
  num_guilds <- 5
  sim <- DAISIE_sim_constant_rate(
    time = island_age,
    M = n_mainland_species,
    pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
    replicates = replicates,
    num_guilds = num_guilds,
    plot_sims = FALSE,
    verbose = FALSE
  )
  expect_true(is.list(sim))
  expect_true(length(sim) == replicates)
})

test_that("Output is silent for nonoceanic_pars[1] != 0 when
          divdepmodel = 'CS'", {
            n_mainland_species <- 1000
            island_age <- 0.4
            clado_rate <- 2.550687345 # cladogenesis rate
            ext_rate <- 2.683454548 # extinction rate
            clade_carr_cap <- 10.0  # clade-level carrying capacity
            imm_rate <- 0.00933207 # immigration rate
            ana_rate <- 1.010073119 # anagenesis rate
            nonoceanic_pars <- c(0.1, 0.9)
            expect_silent(
              DAISIE_sim_constant_rate(
                time = island_age,
                M = n_mainland_species,
                pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
                replicates = 1,
                nonoceanic_pars = nonoceanic_pars,
                plot_sims = FALSE,
                verbose = FALSE
              )
            )
          })

test_that("output is correct for nonoceanic_pars[1] != 0 when
          divdepmodel = 'CS'", {
            n_mainland_species <- 1000
            island_age <- 0.4
            clado_rate <- 2.550687345 # cladogenesis rate
            ext_rate <- 2.683454548 # extinction rate
            clade_carr_cap <- 10.0  # clade-level carrying capacity
            imm_rate <- 0.00933207 # immigration rate
            ana_rate <- 1.010073119 # anagenesis rate
            replicates <- 1
            nonoceanic_pars <- c(0.1, 0.9)
            sim <- DAISIE_sim_constant_rate(time = island_age,
                                            M = n_mainland_species,
                                            pars = c(clado_rate,
                                                     ext_rate,
                                                     clade_carr_cap,
                                                     imm_rate,
                                                     ana_rate),
                                            replicates = replicates,
                                            divdepmodel = "CS",
                                            nonoceanic_pars = nonoceanic_pars,
                                            plot_sims = FALSE,
                                            verbose = FALSE)
            expect_true(is.list(sim))
            expect_true(length(sim) == replicates)
          })

test_that("Output is silent for nonoceanic_pars[1] != 0 when
          divdepmodel = 'IW'", {
            n_mainland_species <- 1000
            island_age <- 0.4
            clado_rate <- 2.550687345 # cladogenesis rate
            ext_rate <- 2.683454548 # extinction rate
            clade_carr_cap <- 10.0  # clade-level carrying capacity
            imm_rate <- 0.00933207 # immigration rate
            ana_rate <- 1.010073119 # anagenesis rate
            nonoceanic_pars <- c(0.1, 0.9)
            expect_silent(
              DAISIE_sim_constant_rate(
                time = island_age,
                M = n_mainland_species,
                pars = c(clado_rate,
                         ext_rate,
                         clade_carr_cap,
                         imm_rate,
                         ana_rate),
                replicates = 1,
                divdepmodel = "IW",
                nonoceanic_pars = nonoceanic_pars,
                plot_sims = FALSE,
                verbose = FALSE
              )
            )
          })

test_that("output is correct for nonoceanic_pars[1] != 0 when
          divdepmodel = 'IW'", {
            n_mainland_species <- 1000
            island_age <- 0.4
            clado_rate <- 2.550687345 # cladogenesis rate
            ext_rate <- 2.683454548 # extinction rate
            clade_carr_cap <- 10.0  # clade-level carrying capacity
            imm_rate <- 0.00933207 # immigration rate
            ana_rate <- 1.010073119 # anagenesis rate
            replicates <- 1
            nonoceanic_pars <- c(0.1, 0.9)
            sim <- DAISIE_sim_constant_rate(time = island_age,
                                            M = n_mainland_species,
                                            pars = c(clado_rate,
                                                     ext_rate,
                                                     clade_carr_cap,
                                                     imm_rate,
                                                     ana_rate),
                                            replicates = replicates,
                                            divdepmodel = "IW",
                                            nonoceanic_pars = nonoceanic_pars,
                                            plot_sims = FALSE,
                                            verbose = FALSE)
            expect_true(is.list(sim))
            expect_true(length(sim) == replicates)
          })

test_that("A non-oceanic run should have native species on the island", {
  n_mainland_species <- 1000
  island_age <- 0.4
  clado_rate <- 2.550687345 # cladogenesis rate
  ext_rate <- 2.683454548 # extinction rate
  clade_carr_cap <- 10.0  # clade-level carrying capacity
  imm_rate <- 0.00933207 # immigration rate
  ana_rate <- 1.010073119 # anagenesis rate
  nonoceanic_pars <- c(0.5, 0.9)
  sim <- DAISIE_sim_constant_rate(
    time = island_age,
    M = n_mainland_species,
    pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
    replicates = 1,
    nonoceanic_pars = nonoceanic_pars,
    plot_sims = FALSE,
    verbose = FALSE
  )
  #number of immigrants (nonendemics) is greater than zero
  expect_gt(sim[[1]][[1]]$stt_all[1, 2], 0)
  #number of anagenetic species (endemic) is greater than zero
  expect_gt(sim[[1]][[1]]$stt_all[1, 3], 0)
})

test_that("Oceanic and non-oceanic should give same results when
          initial sampling is zero", {
            n_mainland_species <- 1000
            island_age <- 0.4
            clado_rate <- 2.550687345 # cladogenesis rate
            ext_rate <- 2.683454548 # extinction rate
            clade_carr_cap <- 10.0  # clade-level carrying capacity
            imm_rate <- 0.00933207 # immigration rate
            ana_rate <- 1.010073119 # anagenesis rate
            set.seed(17)
            oceanic_sim <- DAISIE_sim_constant_rate(
              time = island_age,
              M = n_mainland_species,
              pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
              replicates = 1,
              plot_sims = FALSE,
              verbose = FALSE
            )
            set.seed(17)
            nonoceanic_sim <- DAISIE_sim_constant_rate(
              time = island_age,
              M = n_mainland_species,
              pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
              replicates = 1,
              nonoceanic_pars = c(0.0, 0.9),
              plot_sims = FALSE,
              verbose = FALSE
            )
            expect_true(all(names(oceanic_sim) == names(nonoceanic_sim)))
          })

test_that("abuse: error when mainland n is not multiple of guild number", {
  n_mainland_species <- 1000
  island_age <- 0.4
  clado_rate <- 2.550687345 # cladogenesis rate
  ext_rate <- 2.683454548 # extinction rate
  clade_carr_cap <- 10.0  # clade-level carrying capacity
  imm_rate <- 0.00933207 # immigration rate
  ana_rate <- 1.010073119 # anagenesis rate
  num_guilds <- 33
  expect_error(
    DAISIE_sim_constant_rate(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      divdepmodel = "GW",
      num_guilds = num_guilds,
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
})

test_that("abuse GW parsing errors", {
  expect_error()
  n_mainland_species <- 100
  island_age <- 0.4
  clado_rate <- 1.0
  ext_rate <- 1.0
  clade_carr_cap <- 10.0
  imm_rate <- 0.01
  ana_rate <- 1.0
  num_guilds <- "nonsense"
  expect_error(
    DAISIE_sim_constant_rate(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      divdepmodel = "GW",
      num_guilds = num_guilds,
      plot_sims = FALSE,
      verbose = FALSE
    ),
    regexp = "num_guilds must be numeric"
  )
})

test_that("abuse IW with more than 5 parameters", {
  n_mainland_species <- 100
  island_age <- 0.4
  expect_error(
    DAISIE_sim_constant_rate(
      time = island_age,
      M = n_mainland_species,
      pars = c(1.0, 1.0, 10.0, 0.01, 1.0, 1.0),
      replicates = 1,
      divdepmodel = "IW",
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
})

test_that("constant rate oceanic CS prints correct output when
          verbose == TRUE", {
            totaltime <- 1
            mainland_n <- 1
            pars <- c(0.4, 0.2, 10, 2, 0.8)
            replicates <- 1
            verbose <- TRUE
            set.seed(1)
            expect_output(
              sim <- DAISIE::DAISIE_sim_constant_rate(time = totaltime,
                                                      M = mainland_n,
                                                      pars = pars,
                                                      replicates = replicates,
                                                      plot_sims = FALSE,
                                                      verbose = verbose),
              regexp = "Island replicate 1"
            )
          })

test_that("constant rate oceanic IW prints correct output when
          verbose == TRUE", {
            totaltime <- 1
            mainland_n <- 2
            pars <- c(0.4, 0.2, 10, 2, 0.8)
            replicates <- 1
            verbose <- TRUE
            set.seed(1)
            expect_output(
              sim <- DAISIE::DAISIE_sim_constant_rate(time = totaltime,
                                                      M = mainland_n,
                                                      pars = pars,
                                                      replicates = replicates,
                                                      plot_sims = FALSE,
                                                      verbose = verbose,
                                                      divdepmodel = "IW"),
              regexp = "Island replicate 1"
            )
          })


test_that("2 type simulation with divdepmodel = 'CS' verbose run should
          print", {
            n_mainland_species <- 100
            island_age <- 0.4
            clado_rate_type_1 <- 1.0
            ext_rate_type_1 <- 1.0
            clade_carr_cap_type_1 <- 10.0
            imm_rate_type_1 <- 0.01
            ana_rate_type_1 <- 1.0
            clado_rate_type_2 <- 1.0
            ext_rate_type_2 <- 1.0
            clade_carr_cap_type_2 <- 10.0
            imm_rate_type_2 <- 0.01
            ana_rate_type_2 <- 1.0
            prop_type2_pool <- 0.1
            replicates_apply_type2 <- FALSE
            expect_output(
              sim <- DAISIE_sim_constant_rate(
                time = island_age,
                M = n_mainland_species,
                pars = c(clado_rate_type_1,
                         ext_rate_type_1,
                         clade_carr_cap_type_1,
                         imm_rate_type_1,
                         ana_rate_type_1,
                         clado_rate_type_2,
                         ext_rate_type_2,
                         clade_carr_cap_type_2,
                         imm_rate_type_2,
                         ana_rate_type_2),
                replicates = 1,
                prop_type2_pool = prop_type2_pool,
                replicates_apply_type2 = replicates_apply_type2,
                plot_sims = FALSE,
                verbose = TRUE
              ),
              regexp = "Island replicate 1"
            )
          })

test_that("A divdepmodel = 'GW' run with verbose should print", {
  n_mainland_species <- 100
  island_age <- 0.4
  clado_rate <- 1.0
  ext_rate <- 1.0
  clade_carr_cap <- 10.0
  imm_rate <- 0.01
  ana_rate <- 1.0
  num_guilds <- 5
  expect_output(
    DAISIE_sim_constant_rate(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      divdepmodel = "GW",
      num_guilds = num_guilds,
      plot_sims = FALSE,
      verbose = TRUE
    ),
    regexp = "Island replicate 1"
  )
})

test_that("2 type, no geodynamics, nonoceanic should give error", {
  pars <- c(0.4, 0.1, 10, 1, 0.5, 0.4, 0.1, 10, 1, 0.5)
  totaltime <- 5
  M <- 10
  verbose <- FALSE
  replicates <- 1
  set.seed(1)
  prop_type2_pool <- 0.4
  nonoceanic_pars <- c(0.5, 0.5)
  expect_error(DAISIE::DAISIE_sim_constant_rate(
    time = totaltime,
    M = M,
    pars = pars,
    replicates = replicates,
    prop_type2_pool = prop_type2_pool,
    nonoceanic_pars = nonoceanic_pars,
    verbose = FALSE)
    )
})
