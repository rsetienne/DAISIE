## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
 collapse = TRUE,
 comment = "#>"
)

## ------------------------------------------------------------------------
library(DAISIE)

## ------------------------------------------------------------------------
n_mainland_species <- 1000
island_age <- 4

## ------------------------------------------------------------------------
n_replicates <- 10

## ------------------------------------------------------------------------
set.seed(42)

## ----fig.width=7, fig.height=7-------------------------------------------
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
  verbose = FALSE
)
DAISIE_plot_sims(island_replicates = island_replicates, use_dev_new = FALSE)

## ----fig.width=7, fig.height=7-------------------------------------------
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
DAISIE_plot_sims(island_replicates_K, use_dev_new = FALSE)

## ----fig.width=7, fig.height=7-------------------------------------------
clado_rate_1 <- 0.38 # cladogenesis rate
ext_rate_1 <- 0.55 # extinction rate
clade_carr_cap_1 <- Inf  # clade-level carrying capacity
imm_rate_1 <- 0.04 # immigration rate
ana_rate_1 <- 1.10 # anagenesis rate

clado_rate_2 <- 0.38 # cladogenesis rate
ext_rate_2 <- ext_rate_1 # extinction rate
clade_carr_cap_2 <- clade_carr_cap_1  # clade-level carrying capacity
imm_rate_2 <- imm_rate_1 # immigration rate
ana_rate_2 <- ana_rate_1 # anagenesis rate

island_replicates_2types <- DAISIE_sim( 
  time = island_age,
  M = n_mainland_species, 
  pars = c(
    clado_rate_1, ext_rate_1, clade_carr_cap_1, imm_rate_1, ana_rate_1,
    clado_rate_2, ext_rate_2, clade_carr_cap_2, imm_rate_2, ana_rate_2
  ),
  replicates = n_replicates, 
  prop_type2_pool = 0.163,
  plot_sims = FALSE,
  verbose = FALSE
)
DAISIE_plot_sims(island_replicates_2types, use_dev_new = FALSE)
