context("DAISIE_sim")

test_that("A clean classic run should produce no output", {
  n_mainland_species <- 1000
  island_age <- 4
  clado_rate <- 2.550687345 # cladogenesis rate
  ext_rate <- 2.683454548 # extinction rate
  clade_carr_cap <- 10.0  # clade-level carrying capacity
  imm_rate <- 0.00933207 # immigration rate
  ana_rate <- 1.010073119 # anagenesis rate
  expect_silent(
    DAISIE_sim( 
      time = island_age, 
      M = n_mainland_species, 
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
})

test_that("A clean ontogeny run should produce no output", {
  n_mainland_species <- 1000
  island_age <- 4
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
  island_ontogeny <- "quadratic"
  extcutoff <- 1000
  
  expect_silent(
    DAISIE_sim(
      time = island_age, 
      M = n_mainland_species, 
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1, 
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
  island_age <- 4
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
  island_ontogeny <- "quadratic"
  extcutoff <- 1000
  keep_final_state = TRUE
  
  expect_silent(
    out <- DAISIE_sim(
      time = island_age, 
      M = n_mainland_species, 
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1, 
      island_ontogeny = island_ontogeny,
      Apars = create_area_params(max_area, peak_time, sharpness, total_island_age),
      Epars = c(mu_min, mu_max),
      extcutoff = extcutoff,
      plot_sims = FALSE,
      verbose = FALSE,
      keep_final_state = keep_final_state
    )
  )
  expect_true(is.matrix(out[[1]][[2]]$island_spec))
})


