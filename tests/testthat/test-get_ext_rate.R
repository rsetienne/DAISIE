context("get ext_rate")

test_that("ext rate is a number", {

  
  get_ext_rate(
    timeval = 0,
    mu = 2,
    Apars = create_area_params(10, 0.5, 1, 15),
    Epars = c(1, 10), 
    island_ontogeny = "quadratic", 
    extcutoff = 1000,
    island_spec = c(),
    K = 10
  )
  
  expect_silent(
    is.numeric(
      get_ext_rate(
        timeval = 0,
        mu = 2,
        Apars = create_area_params(10, 0.5, 1, 15),
        Epars = c(1, 10), 
        island_ontogeny = "quadratic", 
        extcutoff = 1000,
        island_spec = c(),
        K = 10
      )
    )
  )
})

test_that("classic behaviour", {
  
  
  carr_cap <- 10
  ps_ext_rate <- 2
  n_species <- 5
  n_mainland_species <- 2
  expected <- DAISIE_calc_clade_ext_rate(
    ps_ext_rate = ps_ext_rate,
    n_species = n_species
  )
  created <- get_ext_rate(
    timeval = 1.0,
    mu = ps_ext_rate,
    Apars =  NULL,
    Epars = NULL,
    island_ontogeny = NULL,
    island_spec = matrix(data = NA, nrow = n_species, ncol = 1),
    K = carr_cap,
    extcutoff = 1000
  )
  expect_equal(expected, created)
})
