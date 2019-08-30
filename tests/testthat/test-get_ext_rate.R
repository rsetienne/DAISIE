context("get ext_rate")

test_that("ext rate is a number", {

  
  get_ext_rate(
    timeval = 0,
    mu = 2,
    Tpars = NULL,
    Apars = create_area_params(10, 0.5, 1, 15),
    Epars = c(1, 10), 
    island_ontogeny = translate_island_ontogeny("beta"),
    extcutoff = 1000,
    island_spec = c(),
    K = 10
  )
  
  expect_silent(
    is.numeric(
      get_ext_rate(
        timeval = 0,
        mu = 2,
        Tpars = NULL,
        Apars = create_area_params(10, 0.5, 1, 15),
        Epars = c(1, 10),
        island_ontogeny = translate_island_ontogeny("const"),
        extcutoff = 1000,
        island_spec = c(),
        K = 10
      )
    )
  )
})
test_that("ext rate is a list", {
  
  Tpars = list(trans_rate = 0.5, 
               immig_rate2 = 0.1, 
               ext_rate2 = 0.2, 
               ana_rate2 = 0.3, 
               clado_rate2 = 0.4, 
               trans_rate2 = 0.5, 
               M2 = 1000)
  get_ext_rate(
    timeval = 0,
    mu = 2,
    Tpars = Tpars,
    Apars = NULL,
    Epars = NULL, 
    island_ontogeny = translate_island_ontogeny("const"),
    extcutoff = 1000,
    island_spec = c(),
    K = 10
  )
  
  expect_silent(
    is.list(
      get_ext_rate(
        timeval = 0,
        mu = 2,
        Tpars = Tpars,
        Apars = NULL,
        Epars = NULL,
        island_ontogeny = translate_island_ontogeny("const"),
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
    Tpars = NULL,
    Apars = NULL,
    Epars = NULL,
    island_ontogeny = translate_island_ontogeny("const"),
    island_spec = matrix(data = NA, nrow = n_species, ncol = 1),
    K = carr_cap,
    extcutoff = 1000
  )
  expect_equal(expected, created)
})
