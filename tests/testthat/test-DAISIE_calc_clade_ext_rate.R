context("DAISIE_calc_clade_ext_rate")

test_that("use", {
  ps_ext_rate <- 1.2
  n_species <- 42
  created <- DAISIE_calc_clade_ext_rate(
    ps_ext_rate = ps_ext_rate,
    n_species = n_species
  )
  expected <- ps_ext_rate * n_species
  expect_equal(created, expected)
})

test_that("abuse", {
  ps_ext_rate <- 1.2
  n_species <- 42
  expect_error(
    DAISIE_calc_clade_ext_rate(
      ps_ext_rate = -123.456,
      n_species = n_species
    )
  )
  expect_error(
    DAISIE_calc_clade_ext_rate(
      ps_ext_rate = ps_ext_rate,
      n_species = -123456
    )
  )
})
