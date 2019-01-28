context("get_immig_rate")

test_that("immig rate plots", {
  immig <- c()
  timepoints <- seq(0, 10, by = 0.01)
  
  for (i in 1:1000) {
    immig[i] <- get_immig_rate(
      timepoints[i], totaltime = 10, gam = 0.001,
       Apars = create_area_params(5000, 0.2, 1, 15), 
       island_spec = matrix(ncol = 1), 
       island_ontogeny = "quadratic", 
       mainland_n = 1000, K = 0.05
    )
  }
  return()

  
})

test_that("classic behavior", {
  
  carr_cap <- 10
  ps_imm_rate <- 0.1
  n_island_species <- 5
  n_mainland_species <- 2
  expected <- DAISIE_calc_clade_imm_rate(
    ps_imm_rate = ps_imm_rate,
    n_island_species = n_island_species,
    n_mainland_species = n_mainland_species,
    carr_cap = carr_cap
  )
  created <- get_immig_rate(
    timeval = 1.0,
    totaltime = 10.0,
    gam = ps_imm_rate,
    Apars =  NULL,
    island_ontogeny = NULL,
    island_spec = matrix(data = NA, nrow = n_island_species, ncol = 1),
    K = carr_cap,
    mainland_n = n_mainland_species
  )
  expect_equal(expected, created)
  
})
