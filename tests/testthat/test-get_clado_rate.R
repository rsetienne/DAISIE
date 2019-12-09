context("get_clado_rate")

test_that("use area constant diversity-independent without hyper_pars", {
  expect_silent(
    get_clado_rate(timeval = 0,
                   lac = 2,
                   ddmodel_sim = 0,
                   hyper_pars = NULL,
                   area_pars = create_area_pars(1, 0, 0, 0, 0, 0),
                   dist_pars = 1,
                   island_ontogeny = translate_island_ontogeny("const"),
                   sea_level = translate_sea_level("const"),
                   num_spec = 0,
                   K = 10
    )
  )
})

test_that("use area constant diversity-independent with hyper_pars", {
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

test_that("use area constant diversity-dependent without hyper_pars", {
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

test_that("use area constant diversity-dependent with hyper_pars", {})

test_that("use area variable (ontogeny) diversity-dependent without
          hyper_pars",{})

test_that("use area variable (ontogeny) diversity-dependent with
          hyper_pars",{})

test_that("use area variable (sea-level) diversity-dependent without
          hyper_pars", {})

test_that("use area variable (sea-level) diversity-dependent with
          hyper_pars", {})

test_that("use area variable (ontogeny and sea-level) diversity-dependent
          without hyper_pars", {})

test_that("use area variable (ontogeny and sea-level) diversity-dependent
          with hyper_pars", {})
