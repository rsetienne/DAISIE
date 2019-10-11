context("DAISIE_sim")

test_that("A clean classic run should produce no output", {
  n_mainland_species <- 1000
  island_age <- 0.4
  clado_rate <- 2.550687345 # cladogenesis rate
  ext_rate <- 2.683454548 # extinction rate
  clade_carr_cap <- 10.0  # clade-level carrying capacity
  imm_rate <- 0.00933207 # immigration rate
  ana_rate <- 1.010073119 # anagenesis rate
  expect_silent(
    DAISIE_sim(
      time = island_age,
      M = n_mainland_species,
      Tpars = NULL,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
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
      mainland_params = NULL,
      divdepmodel = 'CS',
      island_type = island_type,
      nonoceanic = NULL,
      prop_type2_pool = NA,
      replicates_apply_type2 = TRUE,
      sample_freq = 25,
      island_ontogeny = island_ontogeny,
      Apars = create_area_params(max_area, peak_time, sharpness, total_island_age),
      Epars = c(mu_min, mu_max),
      Tpars = NULL,
      keep_final_state = FALSE,
      stored_data = NULL,
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
      mainland_params = NULL,
      divdepmodel = 'CS',
      island_type = "oceanic",
      nonoceanic = NULL,
      prop_type2_pool = NA,
      replicates_apply_type2 = TRUE,
      sample_freq = 25,
      island_ontogeny = island_ontogeny,
      Apars = create_area_params(max_area, peak_time, sharpness, total_island_age),
      Epars = c(mu_min, mu_max),
      Tpars = NULL,
      stored_data = NULL,
      plot_sims = FALSE,
      verbose = FALSE,
      keep_final_state = keep_final_state
    )
  )
  expect_true(is.matrix(out[[1]][[2]]$island_spec) || length(out[[1]][[2]]$branching_times) == 1)
})

test_that("A DAISIE IW simulation that produces empty islands works", {

  set.seed(1)
  time <- 10
  M <- 1000
  pars <- c(1.0e+00, 4.0e-01, 1.5e+01, 1.0e-04, 2.0e-01)
  replicates <- 10
  divdepmodel <-  "IW"

  expect_silent(
    out <- DAISIE::DAISIE_sim(
      time = time,
      M = M,
      pars = pars,
      replicates = 10,
      divdepmodel = "IW",
      verbose = FALSE,
      plot_sims = FALSE
    )
  )

  expect_true(is.matrix(out[[1]][[1]]$stt_all) || length(out[[1]][[2]]$brts_table) == 1)
})



test_that("Vignette demo_sim example 1", {
  n_mainland_species <- 1000
  island_age <- 4
  n_replicates <- 10
  set.seed(42)
  clado_rate <- 2.550687345 # cladogenesis rate
  ext_rate <- 2.683454548 # extinction rate
  clade_carr_cap <- Inf # clade-level carrying capacity
  imm_rate <- 0.00933207 # immigration rate
  ana_rate <- 1.010073119 # anagenesis rate

  island_replicates <- DAISIE_sim(
    time = island_age,
    M = n_mainland_species,
    pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
    replicates = n_replicates,
    plot_sims = FALSE,
    verbose = FALSE,
    Apars = NULL
  )




  clado_rate <- 2.550687345 # cladogenesis rate
  ext_rate <- 2.683454548 # extinction rate
  clade_carr_cap <- 10.0  # clade-level carrying capacity
  imm_rate <- 0.00933207 # immigration rate
  ana_rate <- 1.010073119 # anagenesis rate

  island_replicates_K <- DAISIE_sim(
    time = island_age,
    M = n_mainland_species,
    pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
    replicates = n_replicates,
    plot_sims = FALSE,
    verbose = FALSE
  )

})
