context("DAISIE_sim")

test_that("A clean classic run should produce no output", {
  n_mainland_species <- 1000
  island_age <- 0.4
  clado_rate <- 2.550687345 # cladogenesis rate
  ext_rate <- 2.683454548 # extinction rate
  clade_carr_cap <- 10.0  # clade-level carrying capacity
  imm_rate <- 0.00933207 # immigration rate
  ana_rate <- 1.010073119 # anagenesis rate
  divdep <- c('lac', 'gam')
  island_type <- "oceanic"
  expect_silent(
    DAISIE_sim( 
      time = island_age, 
      M = n_mainland_species, 
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      divdep = divdep,
      island_type = island_type,
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
  divdep <- c('lac', 'gam')
  island_type <- "oceanic"
  max_area <- 1000
  peak_time <- 0.1
  sharpness <- 1
  total_island_age <- 10
  mu_min <- 0.5
  mu_max <- 100
  island_ontogeny <- "beta"
  extcutoff <- 1000
  
  expect_silent(
    DAISIE_sim(
      time = island_age, 
      M = n_mainland_species, 
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1, 
      divdep = divdep,
      island_type = island_type,
      island_ontogeny = island_ontogeny,
      Apars = create_area_params(max_area, peak_time, sharpness, total_island_age),
      Epars = c(mu_min, mu_max),
      extcutoff = extcutoff,
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
})

test_that("A keep last final state ontogeny run should produce no output and store island_spec", {
  n_mainland_species <- 1000
  island_age <- 0.4
  clado_rate <- 0.0001 # cladogenesis rate
  ext_rate <- 2.683454548 # extinction rate (not used)
  clade_carr_cap <- 0.05  # clade-level carrying capacity
  imm_rate <- 0.001 # immigration rate
  ana_rate <- 0.1 # anagenesis rate
  divdep <- c('lac', 'gam')
  island_type <- "oceanic"
  max_area <- 1000
  peak_time <- 0.1
  sharpness <- 1
  total_island_age <- 10
  mu_min <- 0.5
  mu_max <- 100
  island_ontogeny <- "beta"
  extcutoff <- 1000
  keep_final_state = TRUE
  
  expect_silent(
    out <- DAISIE_sim(
      time = island_age, 
      M = n_mainland_species, 
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1, 
      divdep = divdep,
      island_type = island_type,
      island_ontogeny = island_ontogeny,
      Apars = create_area_params(max_area, peak_time, sharpness, total_island_age),
      Epars = c(mu_min, mu_max),
      extcutoff = extcutoff,
      plot_sims = FALSE,
      verbose = FALSE,
      keep_final_state = keep_final_state
    )
  )
  expect_true(is.matrix(out[[1]][[2]]$island_spec) || length(out[[1]][[2]]$branching_times) == 1)
})

test_that("An oceanic run with diversity-dependent mu should produce no output", {
  n_mainland_species <- 1000
  island_age <- 0.4
  clado_rate <- 2.550687345 # cladogenesis rate
  ext_rate <- 2.683454548 # extinction rate
  clade_carr_cap <- 10.0  # clade-level carrying capacity
  imm_rate <- 0.00933207 # immigration rate
  ana_rate <- 1.010073119 # anagenesis rate
  divdep <- c('lac', 'mu', 'gam')
  island_type <- "oceanic"
  expect_silent(
    DAISIE_sim( 
      time = island_age, 
      M = n_mainland_species, 
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      divdep = divdep,
      island_type = island_type,
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
})

test_that("An oceanic run with diversity-independent rates should produce no output", {
  n_mainland_species <- 1000
  island_age <- 0.4
  clado_rate <- 2.550687345 # cladogenesis rate
  ext_rate <- 2.683454548 # extinction rate
  clade_carr_cap <- 10.0  # clade-level carrying capacity
  imm_rate <- 0.00933207 # immigration rate
  ana_rate <- 1.010073119 # anagenesis rate
  divdep <- NULL
  island_type <- "oceanic"
  expect_silent(
    DAISIE_sim( 
      time = island_age, 
      M = n_mainland_species, 
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      divdep = divdep,
      island_type = island_type,
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
})


test_that("A non-oceanic run should produce no output", {
  n_mainland_species <- 1000
  island_age <- 0.4
  clado_rate <- 2.550687345 # cladogenesis rate
  ext_rate <- 2.683454548 # extinction rate
  clade_carr_cap <- 10.0  # clade-level carrying capacity
  imm_rate <- 0.00933207 # immigration rate
  ana_rate <- 1.010073119 # anagenesis rate
  divdep <- NULL
  island_type <- "nonoceanic"
  nonoceanic <- c(0.1, 0.9)
  expect_silent(
    DAISIE_sim( 
      time = island_age, 
      M = n_mainland_species, 
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      divdep = divdep,
      island_type = island_type,
      nonoceanic = nonoceanic,
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
})


test_that("A non-oceanic run should have native species on the island", {
  n_mainland_species <- 1000
  island_age <- 0.4
  clado_rate <- 2.550687345 # cladogenesis rate
  ext_rate <- 2.683454548 # extinction rate
  clade_carr_cap <- 10.0  # clade-level carrying capacity
  imm_rate <- 0.00933207 # immigration rate
  ana_rate <- 1.010073119 # anagenesis rate
  divdep <- NULL
  island_type <- "nonoceanic"
  nonoceanic <- c(0.1, 0.9)
  sim <- DAISIE_sim( 
      time = island_age, 
      M = n_mainland_species, 
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      divdep = divdep,
      island_type = island_type,
      nonoceanic = nonoceanic,
      plot_sims = FALSE,
      verbose = FALSE
    )
  expect_gt(sim$stt_table[,1], 0)
  expect_gt(sim$stt_tabel[,2], 0)
})

