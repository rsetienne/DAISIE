context("DAISIE_sim_time_dependent")

test_that("A clean CS ontogeny run should produce no output", {
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
  sea_level_amplitude <- 0
  sea_level_frequency <- 0
  island_gradient_angle <- 0
  mu_min <- 0.5
  mu_max <- 100
  island_ontogeny <- "beta"
  sea_level <- "const"
  extcutoff <- 1000
  expect_silent(
    DAISIE_sim_time_dependent(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      area_pars = create_area_pars(max_area,
                                   peak_time,
                                   sharpness,
                                   total_island_age,
                                   sea_level_amplitude,
                                   sea_level_frequency,
                                   island_gradient_angle),
      ext_pars = c(mu_min, mu_max),
      extcutoff = extcutoff,
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
})

test_that("A clean IW ontogeny run should produce no output", {
  n_mainland_species <- 1000
  island_age <- 0.4
  clado_rate <- 0.0001 # cladogenesis rate
  ext_rate <- 2.683454548 # extinction rate (not used)
  clade_carr_cap <- 0.05  # clade-level carrying capacity
  imm_rate <- 0.001 # immigration rate
  ana_rate <- 0.1 # anagenesis rate
  divdepmodel <- "IW"
  max_area <- 1000
  peak_time <- 0.1
  sharpness <- 1
  total_island_age <- 10
  sea_level_amplitude <- 0
  sea_level_frequency <- 0
  island_gradient_angle <- 0
  mu_min <- 0.5
  mu_max <- 100
  island_ontogeny <- "beta"
  sea_level <- "const"
  extcutoff <- 1000
  expect_silent(
    DAISIE_sim_time_dependent(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      divdepmodel = divdepmodel,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      area_pars = create_area_pars(max_area,
                                   peak_time,
                                   sharpness,
                                   total_island_age,
                                   sea_level_amplitude,
                                   sea_level_frequency,
                                   island_gradient_angle),
      ext_pars = c(mu_min, mu_max),
      extcutoff = extcutoff,
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
})

test_that("A clean GW ontogeny run should produce no output", {
  n_mainland_species <- 1000
  island_age <- 0.4
  clado_rate <- 0.0001 # cladogenesis rate
  ext_rate <- 2.683454548 # extinction rate (not used)
  clade_carr_cap <- 0.05  # clade-level carrying capacity
  imm_rate <- 0.001 # immigration rate
  ana_rate <- 0.1 # anagenesis rate
  divdepmodel <- "GW"
  num_guilds <- 10
  max_area <- 1000
  peak_time <- 0.1
  sharpness <- 1
  total_island_age <- 10
  sea_level_amplitude <- 0
  sea_level_frequency <- 0
  island_gradient_angle <- 0
  mu_min <- 0.5
  mu_max <- 100
  island_ontogeny <- "beta"
  sea_level <- "const"
  extcutoff <- 1000
  expect_silent(
    DAISIE_sim_time_dependent(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      divdepmodel = divdepmodel,
      num_guilds = num_guilds,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      area_pars = create_area_pars(max_area,
                                   peak_time,
                                   sharpness,
                                   total_island_age,
                                   sea_level_amplitude,
                                   sea_level_frequency,
                                   island_gradient_angle),
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
  max_area <- 1000
  peak_time <- 0
  sharpness <- 0
  total_island_age <- 10
  sea_level_amplitude <- 50
  sea_level_frequency <- 10
  island_gradient_angle <- 45
  mu_min <- 0.5
  mu_max <- 100
  island_ontogeny <- "const"
  sea_level <- "sine"
  extcutoff <- 1000
  expect_silent(
    out <- DAISIE_sim_time_dependent(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      island_ontogeny = island_ontogeny,
      area_pars = create_area_pars(max_area,
                                   peak_time,
                                   sharpness,
                                   total_island_age,
                                   sea_level_amplitude,
                                   sea_level_frequency,
                                   island_gradient_angle),
      ext_pars = c(mu_min, mu_max),
      extcutoff = extcutoff,
      sea_level = sea_level,
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
})
