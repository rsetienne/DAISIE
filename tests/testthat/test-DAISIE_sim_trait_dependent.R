context("DAISIE_sim_trait_dependent")

test_that("A clean CS two trait simulation run should produce no output", {
  n_mainland_species <- 1000
  island_age <- 0.4
  clado_rate <- 0.0001 # cladogenesis rate
  ext_rate <- 2.683454548 # extinction rate (not used)
  clade_carr_cap <- 0.05  # clade-level carrying capacity
  imm_rate <- 0.001 # immigration rate
  ana_rate <- 0.1 # anagenesis rate
  trans_rate <- 0.5 # transition rate
  trait_pars <- create_trait_pars(trans_rate = trans_rate,
                    immig_rate2 = imm_rate / 2,
                    ext_rate2 = ext_rate * 2,
                    ana_rate2 = ana_rate / 2,
                    clado_rate2 = clado_rate / 2,
                    trans_rate2 = trans_rate / 2,
                    M2 = n_mainland_species / 2)
  island_ontogeny <- "const"
  sea_level <- "const"
  extcutoff <- 1000
  expect_silent(
    DAISIE_sim_trait_dependent(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      trait_pars = trait_pars,
      extcutoff = extcutoff,
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
})

test_that("A clean IW two trait simulation run should produce no output", {
  n_mainland_species <- 1000
  island_age <- 0.4
  clado_rate <- 0.0001 # cladogenesis rate
  ext_rate <- 2.683454548 # extinction rate (not used)
  clade_carr_cap <- 0.05  # clade-level carrying capacity
  imm_rate <- 0.001 # immigration rate
  ana_rate <- 0.1 # anagenesis rate
  trans_rate <- 0.5 # transition rate
  trait_pars <- create_trait_pars(trans_rate = trans_rate,
                                  immig_rate2 = imm_rate / 2,
                                  ext_rate2 = ext_rate * 2,
                                  ana_rate2 = ana_rate / 2,
                                  clado_rate2 = clado_rate / 2,
                                  trans_rate2 = trans_rate / 2,
                                  M2 = n_mainland_species / 2)
  divdepmodel <- "IW"
  island_ontogeny <- "const"
  sea_level <- "const"
  extcutoff <- 1000
  expect_silent(
    DAISIE_sim_trait_dependent(
      time = island_age,
      M = n_mainland_species,
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      trait_pars = trait_pars,
      replicates = 1,
      divdepmodel = divdepmodel,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      extcutoff = extcutoff,
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
})
