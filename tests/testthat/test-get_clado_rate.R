context("get_clado_rate")

test_that("classic behaviour", {
  ps_clado_rate <- 0.2
  carr_cap <- 9
  n_species <- 4
  
  created <- get_clado_rate(
    timeval = 5,
    totaltime = 7.0,
    lac = ps_clado_rate,
    Apars = NULL,
    island_ontogeny = NULL,
    island_spec = matrix(NA, nrow = n_species, ncol = 1),
    K = carr_cap
  )
    
  expected <- DAISIE_calc_clade_clado_rate(
    ps_clado_rate = ps_clado_rate,
    n_species = n_species,
    carr_cap = carr_cap
  )
  expect_equal(created, expected)
})
