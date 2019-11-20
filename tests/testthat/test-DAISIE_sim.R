context("DAISIE_sim")

test_that("A divdepmodel = 'CS' run should produce no output", {
  n_mainland_species <- 100
  island_age <- 0.4
  clado_rate <- 1.0
  ext_rate <- 1.0
  clade_carr_cap <- 10.0
  imm_rate <- 0.01
  ana_rate <- 1.0
  ddmodel_sim <- 11
  island_type <- "oceanic"
  expect_silent(
    DAISIE_sim(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      ddmodel_sim = ddmodel_sim,
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
  ddmodel_sim <- 11
  island_type <- "oceanic"
  expect_silent(
    DAISIE_sim(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      divdepmodel = "IW",
      ddmodel_sim = ddmodel_sim,
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
  ddmodel_sim <- 11
  island_type <- "oceanic"
  num_guilds <- 5
  expect_silent(
    DAISIE_sim(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      divdepmodel = "GW",
      ddmodel_sim = ddmodel_sim,
      island_type = island_type,
      num_guilds = num_guilds,
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
  ddmodel_sim <- 11
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
      ddmodel_sim = ddmodel_sim,
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
  island_age <- 0.4
  clado_rate <- 0.0001 # cladogenesis rate
  ext_rate <- 2.683454548 # extinction rate (not used)
  clade_carr_cap <- 0.05  # clade-level carrying capacity
  imm_rate <- 0.001 # immigration rate
  ana_rate <- 0.1 # anagenesis rate
  ddmodel_sim <- 11
  island_type <- "oceanic"
  max_area <- 1000
  peak_time <- 0
  sharpness <- 0
  total_island_age <- 0.4
  sea_level_amplitude <- 50
  sea_level_frequency <- 10
  mu_min <- 0.5
  mu_max <- 100
  island_ontogeny <- "const"
  sea_level <- "sine"
  extcutoff <- 1000
  expect_silent(
    DAISIE_sim(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      ddmodel_sim = ddmodel_sim,
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
    ddmodel_sim <- 11
    island_type <- "oceanic"
    sim <- DAISIE_sim(
        time = island_age,
        M = n_mainland_species,
        pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
        replicates = replicates,
        ddmodel_sim = ddmodel_sim,
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
    ddmodel_sim <- 11
    island_type <- "oceanic"
    sim <- DAISIE_sim(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = replicates,
      divdepmodel = "IW",
      ddmodel_sim = ddmodel_sim,
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
  ddmodel_sim <- 11
  island_type <- "oceanic"
  num_guilds <- 5
  sim <- DAISIE_sim(
    time = island_age,
    M = n_mainland_species,
    pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
    replicates = replicates,
    ddmodel_sim = ddmodel_sim,
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
  ddmodel_sim <- 11
  island_type <- "oceanic"
  expect_silent(
    DAISIE_sim(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      ddmodel_sim = ddmodel_sim,
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
  ddmodel_sim <- 11
  island_type <- "oceanic"
  expect_silent(
    DAISIE_sim(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      ddmodel_sim = ddmodel_sim,
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
  ddmodel_sim <- 11
  island_type <- "nonoceanic"
  nonoceanic_pars <- c(0.1, 0.9)
  expect_silent(
    DAISIE_sim(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      ddmodel_sim = ddmodel_sim,
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
            ddmodel_sim <- 11
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
                              ddmodel_sim = ddmodel_sim,
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
            ddmodel_sim <- 11
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
                ddmodel_sim = ddmodel_sim,
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
            ddmodel_sim <- 11
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
                              ddmodel_sim = ddmodel_sim,
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
  ddmodel_sim <- 11
  island_type <- "nonoceanic"
  nonoceanic_pars <- c(0.5, 0.9)
  sim <- DAISIE_sim(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      ddmodel_sim = ddmodel_sim,
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
  ddmodel_sim <- 11
  island_type <- "oceanic"
  num_guilds <- 33
  expect_error(
    DAISIE_sim(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      divdepmodel = "GW",
      ddmodel_sim = ddmodel_sim,
      island_type = island_type,
      num_guilds = num_guilds,
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
})

test_that("multi-K CS models gives different Ks for each clade", {
  n_mainland_species <- 50
  island_age <- 10.0
  clado_rate <- 2.0
  ext_rate <- 2.0
  clade_carr_cap <- 20.0
  imm_rate <- 0.1
  ana_rate <- 1.0
  ddmodel_sim <- 11
  island_type <- "oceanic"
  k_dist_pars <- c(3, 0.5)
  set.seed(2)
  sim <- DAISIE_sim(
    time = island_age,
    M = n_mainland_species,
    pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
    replicates = 1,
    ddmodel_sim = ddmodel_sim,
    island_type = island_type,
    k_dist_pars = k_dist_pars,
    plot_sims = FALSE,
    verbose = FALSE
    )
  expected_carrying_capacities <- c(2.565935, 11.281834, 8.529117, 3.835293,
                                  12.224073, 2.390876, 2.564162, 4.239230,
                                  3.380080, 6.099572, 3.567437, 3.936415,
                                  11.342538 , 3.279662, 10.964393, 4.810348,
                                  3.492220, 10.256807, 6.489005, 4.900999,
                                  14.528565, 2.801281, 1.427115, 6.076641,
                                  4.599708, 9.270274, 9.412968, 5.319512,
                                  2.918186, 16.762599, 4.719899, 7.534042,
                                  5.824977, 2.181839, 3.538386, 3.306739,
                                  5.097120, 4.932397, 7.536610, 18.117815,
                                  7.593945, 20.128071, 2.994530, 10.467172,
                                  5.689259, 5.654756, 8.055567, 9.116205,
                                  5.114849, 7.788854)
  simulated_carrying_capacities <- sim[[1]][[5]]$all_carrying_capacities
  expect_true(all.equal(expected_carrying_capacities,
                        simulated_carrying_capacities,
                        tolerance = 1e-7))
})

test_that("A multi-K GW models gives different Ks for each clade", {
  n_mainland_species <- 50
  island_age <- 10.0
  clado_rate <- 2.0
  ext_rate <- 2.0
  clade_carr_cap <- 20.0
  imm_rate <- 0.1
  ana_rate <- 1.0
  ddmodel_sim <- 11
  island_type <- "oceanic"
  num_guilds <- 5
  k_dist_pars <- c(3, 0.5)
  set.seed(2)
  sim <- DAISIE_sim(
    time = island_age,
    M = n_mainland_species,
    pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
    replicates = 1,
    divdepmodel = "GW",
    ddmodel_sim = ddmodel_sim,
    island_type = island_type,
    num_guilds = num_guilds,
    k_dist_pars = k_dist_pars,
    plot_sims = FALSE,
    verbose = FALSE
    )
  expected_carrying_capacities <- c(2.565935, 4.239230, 2.278893, 6.339452, 5.824977)
  simulated_carrying_capacities <- sim[[1]][[7]]$all_carrying_capacities
  expect_true(all.equal(expected_carrying_capacities,
                        simulated_carrying_capacities,
                        tolerance = 1e-7))
  })
