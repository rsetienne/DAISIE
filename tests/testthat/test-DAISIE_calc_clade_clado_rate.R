context("DAISIE_calc_clade_clado_rate")

test_that("use, below carrying capacity", {
  ps_clado_rate <- 12.34
  n_species <- 42
  carr_cap <- 314
  created <- DAISIE_calc_clade_clado_rate(
    ps_clado_rate = ps_clado_rate,
    n_species = n_species,
    carr_cap = carr_cap
  )
  expected <- n_species * ps_clado_rate * (1.0 - (n_species / carr_cap))
  expect_equal(created, expected)
})

test_that("use, above carrying capacity", {
  ps_clado_rate <- 12.34
  carr_cap <- 314
  n_species <- carr_cap + 42
  created <- DAISIE_calc_clade_clado_rate(
    ps_clado_rate = ps_clado_rate,
    n_species = n_species,
    carr_cap = carr_cap
  )
  expected <- 0.0
  expect_equal(created, expected)
})

test_that("abuse", {
  ps_clado_rate <- 12.34
  n_species <- 42
  carr_cap <- 314
  expect_error(
    DAISIE_calc_clade_clado_rate(
      ps_clado_rate = -12.34,
      n_species = n_species,
      carr_cap = carr_cap
    )
  )
  expect_error(
    DAISIE_calc_clade_clado_rate(
      ps_clado_rate = ps_clado_rate,
      n_species = -1234,
      carr_cap = carr_cap
    )
  )
  expect_error(
    DAISIE_calc_clade_clado_rate(
      ps_clado_rate = ps_clado_rate,
      n_species = n_species,
      carr_cap = -1234
    )
  )
})
