test_that("A clean CS ontogeny run should produce no output", {
  n_mainland_species <- 1000
  island_age <- 0.4
  clado_rate <- 0.0001 # cladogenesis rate
  ext_rate <- 0.5 # extinction rate
  clade_carr_cap <- 0.05  # clade-level carrying capacity
  imm_rate <- 0.001 # immigration rate
  ana_rate <- 0.1 # anagenesis rate
  max_area <- 1000
  current_area <- 500
  peak_time <- 0.1
  sharpness <- 1
  total_island_age <- 10
  sea_level_amplitude <- 0
  sea_level_frequency <- 0
  island_gradient_angle <- 0
  island_ontogeny <- "beta"
  sea_level <- "const"
  extcutoff <- 1000
  hyper_pars <- create_hyper_pars(d = 0.2, x = 0.1)
  expect_silent(
    DAISIE_sim_time_dep(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      area_pars = create_area_pars(
        max_area = max_area,
        current_area = current_area,
        proportional_peak_t = peak_time,
        total_island_age = total_island_age,
        sea_level_amplitude = sea_level_amplitude,
        sea_level_frequency = sea_level_frequency,
        island_gradient_angle = island_gradient_angle),
      hyper_pars = hyper_pars,
      extcutoff = extcutoff,
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
})



test_that("use CS split-rates with cond", {
  set.seed(Sys.time()) # Always run a different sim
  n_mainland_species <- 1000
  island_age <- 0.4
  clado_rate <- 0.001 # cladogenesis rate
  ext_rate <- 0.5 # extinction rate
  clade_carr_cap <- 0.05  # clade-level carrying capacity
  imm_rate <- 0.01 # immigration rate
  ana_rate <- 0.1 # anagenesis rate
  max_area <- 1000
  current_area <- 500
  peak_time <- 0.1
  sharpness <- 1
  total_island_age <- 10
  sea_level_amplitude <- 0
  sea_level_frequency <- 0
  island_gradient_angle <- 0
  island_ontogeny <- "beta"
  sea_level <- "const"
  extcutoff <- 1000
  hyper_pars <- create_hyper_pars(d = 0.2, x = 0.1)
  cond <- 5
  expect_silent(
    out <- DAISIE_sim_time_dep(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      area_pars = create_area_pars(
        max_area = max_area,
        current_area = current_area,
        proportional_peak_t = peak_time,
        total_island_age = total_island_age,
        sea_level_amplitude = sea_level_amplitude,
        sea_level_frequency = sea_level_frequency,
        island_gradient_angle = island_gradient_angle),
      hyper_pars = hyper_pars,
      extcutoff = extcutoff,
      plot_sims = FALSE,
      verbose = FALSE,
      cond = cond
    )
  )

  expect_true(out[[1]][[1]]$stt_all[nrow(out[[1]][[1]]$stt_all), 5] >= cond)
})


test_that("CS split-rates with cond to without", {
  set.seed(1) # Always run a different sim
  n_mainland_species <- 1000
  island_age <- 0.4
  clado_rate <- 0.001 # cladogenesis rate
  ext_rate <- 0.5 # extinction rate
  clade_carr_cap <- 0.05  # clade-level carrying capacity
  imm_rate <- 0.01 # immigration rate
  ana_rate <- 0.1 # anagenesis rate
  max_area <- 1000
  current_area <- 500
  peak_time <- 0.1
  sharpness <- 1
  total_island_age <- 10
  sea_level_amplitude <- 0
  sea_level_frequency <- 0
  island_gradient_angle <- 0
  island_ontogeny <- "beta"
  sea_level <- "const"
  extcutoff <- 1000
  hyper_pars <- create_hyper_pars(d = 0.2, x = 0.1)
  cond <- 5
  expect_silent(
    out_cond <- DAISIE_sim_time_dep(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      area_pars = create_area_pars(
        max_area = max_area,
        current_area = current_area,
        proportional_peak_t = peak_time,
        total_island_age = total_island_age,
        sea_level_amplitude = sea_level_amplitude,
        sea_level_frequency = sea_level_frequency,
        island_gradient_angle = island_gradient_angle),
      hyper_pars = hyper_pars,
      extcutoff = extcutoff,
      plot_sims = FALSE,
      verbose = FALSE,
      cond = cond
    )
  )
  expect_true(
    out_cond[[1]][[1]]$stt_all[nrow(out_cond[[1]][[1]]$stt_all), 5] >= cond
  )

  set.seed(1)
  expect_silent(
    out_no_cond <- DAISIE_sim_time_dep(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      area_pars = create_area_pars(
        max_area = max_area,
        current_area = current_area,
        proportional_peak_t = peak_time,
        total_island_age = total_island_age,
        sea_level_amplitude = sea_level_amplitude,
        sea_level_frequency = sea_level_frequency,
        island_gradient_angle = island_gradient_angle),
      hyper_pars = hyper_pars,
      extcutoff = extcutoff,
      plot_sims = FALSE,
      verbose = FALSE,
      cond = 0
    )
  )
  expect_true(
    out_no_cond[[1]][[1]]$stt_all[nrow(out_no_cond[[1]][[1]]$stt_all), 5] < cond
    )
})




test_that("A clean IW ontogeny run should produce no output", {
  n_mainland_species <- 1000
  island_age <- 0.4
  clado_rate <- 0.0001 # cladogenesis rate
  ext_rate <- 0.5 # extinction rate
  clade_carr_cap <- 0.05  # clade-level carrying capacity
  imm_rate <- 0.001 # immigration rate
  ana_rate <- 0.1 # anagenesis rate
  divdepmodel <- "IW"
  max_area <- 1000
  current_area <- 500
  peak_time <- 0.1
  total_island_age <- 10
  sea_level_amplitude <- 0
  sea_level_frequency <- 0
  island_gradient_angle <- 0
  island_ontogeny <- "beta"
  sea_level <- "const"
  extcutoff <- 1000
  hyper_pars <- create_hyper_pars(d = 0.2, x = 0.1)
  nonoceanic_pars <- c(0, 0)
  expect_silent(
    DAISIE_sim_time_dep(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      divdepmodel = divdepmodel,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      area_pars = create_area_pars(
        max_area = max_area,
        current_area = current_area,
        proportional_peak_t = peak_time,
        total_island_age = total_island_age,
        sea_level_amplitude = sea_level_amplitude,
        sea_level_frequency = sea_level_frequency,
        island_gradient_angle = island_gradient_angle),
      hyper_pars = hyper_pars,
      nonoceanic_pars = nonoceanic_pars,
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
  ext_rate <- 0.5 # extinction rate
  clade_carr_cap <- 0.05  # clade-level carrying capacity
  imm_rate <- 0.001 # immigration rate
  ana_rate <- 0.1 # anagenesis rate
  divdepmodel <- "GW"
  num_guilds <- 10
  max_area <- 1000
  current_area <- 500
  peak_time <- 0.1
  total_island_age <- 10
  sea_level_amplitude <- 0
  sea_level_frequency <- 0
  island_gradient_angle <- 0
  island_ontogeny <- "beta"
  sea_level <- "const"
  extcutoff <- 1000
  hyper_pars <- create_hyper_pars(d = 0.2, x = 0.1)
  nonoceanic_pars <- c(0, 0)
  expect_silent(
    DAISIE_sim_time_dep(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      divdepmodel = divdepmodel,
      num_guilds = num_guilds,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      area_pars = create_area_pars(
        max_area =  max_area,
        current_area = current_area,
        proportional_peak_t = peak_time,
        total_island_age = total_island_age,
        sea_level_amplitude = sea_level_amplitude,
        sea_level_frequency = sea_level_frequency,
        island_gradient_angle = island_gradient_angle),
      hyper_pars = hyper_pars,
      nonoceanic_pars = nonoceanic_pars,
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
  ext_rate <- 0.5 # extinction rate
  clade_carr_cap <- 0.05  # clade-level carrying capacity
  imm_rate <- 0.001 # immigration rate
  ana_rate <- 0.1 # anagenesis rate
  max_area <- 1000
  current_area <- 500
  peak_time <- 0
  total_island_age <- 10
  sea_level_amplitude <- 50
  sea_level_frequency <- 10
  island_gradient_angle <- 45
  island_ontogeny <- "const"
  sea_level <- "sine"
  extcutoff <- 1000
  hyper_pars <- create_hyper_pars(d = 0.2, x = 0.1)
  nonoceanic_pars <- c(0, 0)
  expect_silent(
    out <- DAISIE_sim_time_dep(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      island_ontogeny = island_ontogeny,
      area_pars = create_area_pars(
        max_area = max_area,
        current_area = current_area,
        proportional_peak_t = peak_time,
        total_island_age = total_island_age,
        sea_level_amplitude = sea_level_amplitude,
        sea_level_frequency = sea_level_frequency,
        island_gradient_angle = island_gradient_angle),
      hyper_pars = hyper_pars,
      nonoceanic_pars = nonoceanic_pars,
      extcutoff = extcutoff,
      sea_level = sea_level,
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
})





test_that("A clean CS ontogeny run with verbose should print to console", {
  n_mainland_species <- 1000
  island_age <- 0.4
  clado_rate <- 0.0001 # cladogenesis rate
  ext_rate <- 0.5 # extinction rate
  clade_carr_cap <- 0.05  # clade-level carrying capacity
  imm_rate <- 0.001 # immigration rate
  ana_rate <- 0.1 # anagenesis rate
  max_area <- 1000
  current_area <- 500
  peak_time <- 0.1
  sharpness <- 1
  total_island_age <- 10
  sea_level_amplitude <- 0
  sea_level_frequency <- 0
  island_gradient_angle <- 0
  island_ontogeny <- "beta"
  sea_level <- "const"
  extcutoff <- 1000
  hyper_pars <- create_hyper_pars(d = 0.2, x = 0.1)
  expect_output(
    DAISIE_sim_time_dep(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      area_pars = create_area_pars(
        max_area = max_area,
        current_area = current_area,
        proportional_peak_t = peak_time,
        total_island_age = total_island_age,
        sea_level_amplitude = sea_level_amplitude,
        sea_level_frequency = sea_level_frequency,
        island_gradient_angle = island_gradient_angle),
      hyper_pars = hyper_pars,
      extcutoff = extcutoff,
      plot_sims = FALSE,
      verbose = TRUE
    ),
    "Island replicate 1"
  )
})

test_that("A clean IW ontogeny run should produce no output", {
  n_mainland_species <- 1000
  island_age <- 0.4
  clado_rate <- 0.0001 # cladogenesis rate
  ext_rate <- 0.5 # extinction rate
  clade_carr_cap <- 0.05  # clade-level carrying capacity
  imm_rate <- 0.001 # immigration rate
  ana_rate <- 0.1 # anagenesis rate
  divdepmodel <- "IW"
  max_area <- 1000
  current_area <- 500
  peak_time <- 0.1
  total_island_age <- 10
  sea_level_amplitude <- 0
  sea_level_frequency <- 0
  island_gradient_angle <- 0
  island_ontogeny <- "beta"
  sea_level <- "const"
  extcutoff <- 1000
  hyper_pars <- create_hyper_pars(d = 0.2, x = 0.1)
  nonoceanic_pars <- c(0, 0)
  expect_output(
    DAISIE_sim_time_dep(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      divdepmodel = divdepmodel,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      area_pars = create_area_pars(
        max_area = max_area,
        current_area = current_area,
        proportional_peak_t = peak_time,
        total_island_age = total_island_age,
        sea_level_amplitude = sea_level_amplitude,
        sea_level_frequency = sea_level_frequency,
        island_gradient_angle = island_gradient_angle),
      hyper_pars = hyper_pars,
      nonoceanic_pars = nonoceanic_pars,
      extcutoff = extcutoff,
      plot_sims = FALSE,
      verbose = TRUE
    ),
    "Island replicate 1"
  )
})

test_that("A clean GW ontogeny run should produce no output", {
  n_mainland_species <- 1000
  island_age <- 0.4
  clado_rate <- 0.0001 # cladogenesis rate
  ext_rate <- 0.5 # extinction rate
  clade_carr_cap <- 0.05  # clade-level carrying capacity
  imm_rate <- 0.001 # immigration rate
  ana_rate <- 0.1 # anagenesis rate
  divdepmodel <- "GW"
  num_guilds <- 10
  max_area <- 1000
  current_area <- 500
  peak_time <- 0.1
  total_island_age <- 10
  sea_level_amplitude <- 0
  sea_level_frequency <- 0
  island_gradient_angle <- 0
  island_ontogeny <- "beta"
  sea_level <- "const"
  extcutoff <- 1000
  hyper_pars <- create_hyper_pars(d = 0.2, x = 0.1)
  nonoceanic_pars <- c(0, 0)
  expect_output(
    DAISIE_sim_time_dep(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      divdepmodel = divdepmodel,
      num_guilds = num_guilds,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      area_pars = create_area_pars(
        max_area =  max_area,
        current_area = current_area,
        proportional_peak_t = peak_time,
        total_island_age = total_island_age,
        sea_level_amplitude = sea_level_amplitude,
        sea_level_frequency = sea_level_frequency,
        island_gradient_angle = island_gradient_angle),
      hyper_pars = hyper_pars,
      nonoceanic_pars = nonoceanic_pars,
      extcutoff = extcutoff,
      plot_sims = FALSE,
      verbose = TRUE
    ),
    "Island replicate 1"
  )
})
test_that("A clean GW ontogeny run should produce no output", {
  n_mainland_species <- 1000
  island_age <- 0.4
  clado_rate <- 0.0001 # cladogenesis rate
  ext_rate <- 0.5 # extinction rate
  clade_carr_cap <- 0.05  # clade-level carrying capacity
  imm_rate <- 0.001 # immigration rate
  ana_rate <- 0.1 # anagenesis rate
  divdepmodel <- "GW"
  num_guilds <- "10"
  max_area <- 1000
  current_area <- 500
  peak_time <- 0.1
  total_island_age <- 10
  sea_level_amplitude <- 0
  sea_level_frequency <- 0
  island_gradient_angle <- 0
  island_ontogeny <- "beta"
  sea_level <- "const"
  extcutoff <- 1000
  hyper_pars <- create_hyper_pars(d = 0.2, x = 0.1)
  nonoceanic_pars <- c(0, 0)
  expect_error(
    DAISIE_sim_time_dep(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      divdepmodel = divdepmodel,
      num_guilds = num_guilds,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      area_pars = create_area_pars(
        max_area =  max_area,
        current_area = current_area,
        proportional_peak_t = peak_time,
        total_island_age = total_island_age,
        sea_level_amplitude = sea_level_amplitude,
        sea_level_frequency = sea_level_frequency,
        island_gradient_angle = island_gradient_angle),
      hyper_pars = hyper_pars,
      nonoceanic_pars = nonoceanic_pars,
      extcutoff = extcutoff,
      plot_sims = FALSE,
      verbose = TRUE
    ),
    "num_guilds must be numeric"
  )
})
