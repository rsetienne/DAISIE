context("get_ext_rate")

test_that("ext rate is a number", {
  expect_silent(
    is.numeric(
      get_ext_rate(
        timeval = 0,
        mu = 2,
        ddmodel_sim = 11,
        hyper_pars = NULL,
        area_pars = create_area_pars(10, 0.5, 1, 15, 0, 0),
        ext_pars = c(1, 10),
        island_ontogeny = translate_island_ontogeny("const"),
        sea_level = translate_sea_level("const"),
        extcutoff = 1000,
        num_spec = 0,
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
    ddmodel_sim = 11,
    hyper_pars = NULL,
    area_pars =  NULL,
    ext_pars = NULL,
    island_ontogeny = translate_island_ontogeny("const"),
    sea_level = translate_sea_level("const"),
    extcutoff = 1000,
    num_spec = n_species,
    K = carr_cap
  )
  expect_equal(expected, created)
})

test_that("use ontogeny behaviour", {
  expect_silent(
    is.numeric(
      DAISIE:::get_ext_rate(
        timeval = 5,
        mu = 2,
        ddmodel_sim = 11,
        hyper_pars = NULL,
        area_pars = create_area_pars(1000, 0.5, 1, 15, 0, 0),
        ext_pars = c(1, 10),
        island_ontogeny = translate_island_ontogeny("beta"),
        sea_level = translate_sea_level("const"),
        extcutoff = 1000,
        num_spec = 10,
        K = 20
      )
    )
  )
})

test_that("use ontogeny behaviour", {
  expect_silent(
    is.numeric(
      DAISIE:::get_ext_rate(
        timeval = 5,
        mu = 2,
        ddmodel_sim = 11,
        hyper_pars = NULL,
        area_pars = create_area_pars(1000, 0.5, 1, 15, 0, 0),
        ext_pars = c(1, 10),
        island_ontogeny = translate_island_ontogeny("beta"),
        sea_level = translate_sea_level("const"),
        extcutoff = 1000,
        num_spec = 10,
        K = 20
      )
    )
  )
})

test_that("use hyper_pars", {
  expect_silent(
    is.numeric(
      DAISIE:::get_ext_rate(
        timeval = 5,
        mu = 2,
        ddmodel_sim = 11,
        hyper_pars = c(1,1,1,1),
        area_pars = create_area_pars(1000, 0.5, 1, 15, 0, 0),
        ext_pars = c(1, 10),
        island_ontogeny = translate_island_ontogeny("beta"),
        sea_level = translate_sea_level("const"),
        extcutoff = 1000,
        num_spec = 10,
        K = 20
      )
    )
  )
})
