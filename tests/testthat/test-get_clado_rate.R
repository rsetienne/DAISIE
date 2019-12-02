context("get_clado_rate")

test_that("classic behaviour", {
  ps_clado_rate <- 0.2
  carr_cap <- 9
  n_species <- 4
  created <- get_clado_rate(
    timeval = 5,
    lac = ps_clado_rate,
    ddmodel_sim = 11,
    hyper_pars = NULL,
    area_pars = NULL,
    dist_pars = NULL,
    island_ontogeny = 0,
    sea_level = 0,
    num_spec = n_species,
    K = carr_cap
  )
  expected <- DAISIE_calc_clade_clado_rate(
    ps_clado_rate = ps_clado_rate,
    n_species = n_species,
    carr_cap = carr_cap
  )
  expect_equal(created, expected)
})
