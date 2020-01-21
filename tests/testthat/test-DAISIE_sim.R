context("DAISIE_sim")

test_that("A divdepmodel = 'CS' run should produce no output", {
  n_mainland_species <- 100
  island_age <- 5
  clado_rate <- 1.0
  ext_rate <- 1.0
  clade_carr_cap <- 10.0
  imm_rate <- 0.01
  ana_rate <- 1.0
  island_type <- "oceanic"
  expect_silent(
    DAISIE_sim(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      island_type = island_type,
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
  island_type <- "oceanic"
  expect_silent(
    DAISIE_sim(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      divdepmodel = "IW",
      island_type = island_type,
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
  island_type <- "oceanic"
  num_guilds <- 5
  expect_silent(
    DAISIE_sim(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      divdepmodel = "GW",
      island_type = island_type,
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
  island_type <- "oceanic"
  prop_type2_pool <- 0.1
  replicates_apply_type2 <- TRUE
  expect_silent(
    sim <- DAISIE_sim(
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
      island_type = island_type,
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
  island_type <- "oceanic"
  prop_type2_pool <- 0.1
  replicates_apply_type2 <- FALSE
  expect_silent(
    sim <- DAISIE_sim(
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
      island_type = island_type,
      prop_type2_pool = prop_type2_pool,
      replicates_apply_type2 = replicates_apply_type2,
      plot_sims = FALSE,
      verbose = FALSE
      )
    )
})

test_that("A clean ontogeny run should produce no output", {
  n_mainland_species <- 1000
  island_age <- 0.4
  clado_rate <- 0.0001 # cladogenesis rate
  ext_rate <- 2.683454548 # extinction rate (not used)
  clade_carr_cap <- 0.05  # clade-level carrying capacity
  imm_rate <- 0.001 # immigration rate
  ana_rate <- 0.1 # anagenesis rate
  island_type <- "oceanic"
  max_area <- 1000
  peak_time <- 0.1
  sharpness <- 1
  total_island_age <- 10
  sea_level_amplitude <- 0
  sea_level_frequency <- 0
  mu_min <- 0.5
  mu_max <- 100
  island_ontogeny <- "beta"
  sea_level <- "const"
  extcutoff <- 1000
  expect_silent(
    DAISIE_sim(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      island_type = island_type,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      area_pars = create_area_pars(max_area,
                                 peak_time,
                                 sharpness,
                                 total_island_age,
                                 sea_level_amplitude,
                                 sea_level_frequency),
      ext_pars = c(mu_min, mu_max),
      extcutoff = extcutoff,
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
})

test_that("A clean sea_level run should produce no output", {
  n_mainland_species <- 1000
  island_age <- 5
  clado_rate <- 0.0001 # cladogenesis rate
  ext_rate <- 2.683454548 # extinction rate (not used)
  clade_carr_cap <- 0.05  # clade-level carrying capacity
  imm_rate <- 0.001 # immigration rate
  ana_rate <- 0.1 # anagenesis rate
  island_type <- "oceanic"
  max_area <- 1000
  peak_time <- 0
  sharpness <- 0
  total_island_age <- 10
  sea_level_amplitude <- 50
  sea_level_frequency <- 10
  mu_min <- 0.5
  mu_max <- 100
  island_ontogeny <- "const"
  sea_level <- "sine"
  extcutoff <- 1000
  expect_silent(
    out <- DAISIE_sim(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      island_type = island_type,
      island_ontogeny = island_ontogeny,
      area_pars = create_area_pars(max_area,
                                   peak_time,
                                   sharpness,
                                   total_island_age,
                                   sea_level_amplitude,
                                   sea_level_frequency),
      ext_pars = c(mu_min, mu_max),
      extcutoff = extcutoff,
      sea_level = sea_level,
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
    island_type <- "oceanic"
    sim <- DAISIE_sim(
        time = island_age,
        M = n_mainland_species,
        pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
        replicates = replicates,
        island_type = island_type,
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
    island_type <- "oceanic"
    sim <- DAISIE_sim(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = replicates,
      divdepmodel = "IW",
      island_type = island_type,
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
  island_type <- "oceanic"
  num_guilds <- 5
  sim <- DAISIE_sim(
    time = island_age,
    M = n_mainland_species,
    pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
    replicates = replicates,
    island_type = island_type,
    num_guilds = num_guilds,
    plot_sims = FALSE,
    verbose = FALSE
  )
  expect_true(is.list(sim))
  expect_true(length(sim) == replicates)
})

test_that("An oceanic run with diversity-dependent mu should produce
          no output", {
  n_mainland_species <- 1000
  island_age <- 0.4
  clado_rate <- 2.550687345 # cladogenesis rate
  ext_rate <- 2.683454548 # extinction rate
  clade_carr_cap <- 10.0  # clade-level carrying capacity
  imm_rate <- 0.00933207 # immigration rate
  ana_rate <- 1.010073119 # anagenesis rate
  island_type <- "oceanic"
  expect_silent(
    DAISIE_sim(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      island_type = island_type,
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
})

test_that("An oceanic run with diversity-independent rates should
          produce no output", {
  n_mainland_species <- 1000
  island_age <- 0.4
  clado_rate <- 2.550687345 # cladogenesis rate
  ext_rate <- 2.683454548 # extinction rate
  clade_carr_cap <- 10.0  # clade-level carrying capacity
  imm_rate <- 0.00933207 # immigration rate
  ana_rate <- 1.010073119 # anagenesis rate
  island_type <- "oceanic"
  expect_silent(
    DAISIE_sim(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      island_type = island_type,
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
})


test_that("Output is silent for island_type = 'nonoceanic' when
          divdepmodel = 'CS'", {
  n_mainland_species <- 1000
  island_age <- 0.4
  clado_rate <- 2.550687345 # cladogenesis rate
  ext_rate <- 2.683454548 # extinction rate
  clade_carr_cap <- 10.0  # clade-level carrying capacity
  imm_rate <- 0.00933207 # immigration rate
  ana_rate <- 1.010073119 # anagenesis rate
  island_type <- "nonoceanic"
  nonoceanic_pars <- c(0.1, 0.9)
  expect_silent(
    DAISIE_sim(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      island_type = island_type,
      nonoceanic_pars = nonoceanic_pars,
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
})

test_that("output is correct for island_type = 'nonoceanic' when
          divdepmodel = 'CS'", {
            n_mainland_species <- 1000
            island_age <- 0.4
            clado_rate <- 2.550687345 # cladogenesis rate
            ext_rate <- 2.683454548 # extinction rate
            clade_carr_cap <- 10.0  # clade-level carrying capacity
            imm_rate <- 0.00933207 # immigration rate
            ana_rate <- 1.010073119 # anagenesis rate
            replicates <- 1
            island_type <- "nonoceanic"
            nonoceanic_pars <- c(0.1, 0.9)
            sim <- DAISIE_sim(time = island_age,
                              M = n_mainland_species,
                              pars = c(clado_rate,
                                       ext_rate,
                                       clade_carr_cap,
                                       imm_rate,
                                       ana_rate),
                              replicates = replicates,
                              divdepmodel = "CS",
                              island_type = island_type,
                              nonoceanic_pars = nonoceanic_pars,
                              plot_sims = FALSE,
                              verbose = FALSE)
            expect_true(is.list(sim))
            expect_true(length(sim) == replicates)
          })

test_that("Output is silent for island_type = 'nonoceanic' when
          divdepmodel = 'IW'", {
            n_mainland_species <- 1000
            island_age <- 0.4
            clado_rate <- 2.550687345 # cladogenesis rate
            ext_rate <- 2.683454548 # extinction rate
            clade_carr_cap <- 10.0  # clade-level carrying capacity
            imm_rate <- 0.00933207 # immigration rate
            ana_rate <- 1.010073119 # anagenesis rate
            island_type <- "nonoceanic"
            nonoceanic_pars <- c(0.1, 0.9)
            expect_silent(
              DAISIE_sim(
                time = island_age,
                M = n_mainland_species,
                pars = c(clado_rate,
                         ext_rate,
                         clade_carr_cap,
                         imm_rate,
                         ana_rate),
                replicates = 1,
                divdepmodel = "IW",
                island_type = island_type,
                nonoceanic_pars = nonoceanic_pars,
                plot_sims = FALSE,
                verbose = FALSE
              )
            )
          })


test_that("output is correct for island_type = 'nonoceanic' when
          divdepmodel = 'IW'", {
            n_mainland_species <- 1000
            island_age <- 0.4
            clado_rate <- 2.550687345 # cladogenesis rate
            ext_rate <- 2.683454548 # extinction rate
            clade_carr_cap <- 10.0  # clade-level carrying capacity
            imm_rate <- 0.00933207 # immigration rate
            ana_rate <- 1.010073119 # anagenesis rate
            replicates <- 1
            island_type <- "nonoceanic"
            nonoceanic_pars <- c(0.1, 0.9)
            sim <- DAISIE_sim(time = island_age,
                              M = n_mainland_species,
                              pars = c(clado_rate,
                                       ext_rate,
                                       clade_carr_cap,
                                       imm_rate,
                                       ana_rate),
                              replicates = replicates,
                              divdepmodel = "IW",
                              island_type = island_type,
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
  island_type <- "nonoceanic"
  nonoceanic_pars <- c(0.5, 0.9)
  sim <- DAISIE_sim(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      island_type = island_type,
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
  oceanic_sim <- DAISIE_sim(
    time = island_age,
    M = n_mainland_species,
    pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
    replicates = 1,
    plot_sims = FALSE,
    verbose = FALSE
  )
  set.seed(17)
  nonoceanic_sim <- DAISIE_sim(
    time = island_age,
    M = n_mainland_species,
    pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
    replicates = 1,
    island_type = "nonoceanic",
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
  island_type <- "oceanic"
  num_guilds <- 33
  expect_error(
    DAISIE_sim(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      divdepmodel = "GW",
      island_type = island_type,
      num_guilds = num_guilds,
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
})

test_that("use split-rates model", {
  expect_silent(
    DAISIE_sim(
      time = 10,
      M = 10,
      pars = c(1, 1, 1, 0.1, 1, 1, 1, 1, 0.1, 1),
      replicates = 1,
      shift_times = 5,
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
})

test_that("abuse split-rates model", {
  expect_error(DAISIE_sim(
    time = 1,
    M = 1,
    pars = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    replicates = 1,
    shift_times = 5,
    verbose = FALSE,
    plot_sims = FALSE
  ))
  expect_error(DAISIE_sim(
    time = 10,
    M = 1,
    pars = c(1, 1, 1, 1, 1),
    replicates = 1,
    shift_times = 5,
    verbose = FALSE,
    plot_sims = FALSE
  ))
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
  island_type <- "oceanic"
  num_guilds <- "nonsense"
  expect_error(
    DAISIE_sim(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      divdepmodel = "GW",
      island_type = island_type,
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
  island_type <- "oceanic"
  expect_error(
    DAISIE_sim(
      time = island_age,
      M = n_mainland_species,
      pars = c(1.0, 1.0, 10.0, 0.01, 1.0, 1.0),
      replicates = 1,
      divdepmodel = "IW",
      island_type = island_type,
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
})
