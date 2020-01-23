context("DAISIE_create_sim_pars")

test_that("output is silent", {
  expect_silent(DAISIE_create_sim_pars())
})

test_that("output is a list", {
  pars <- DAISIE_create_sim_pars()
  expect_true(class(pars) == "list")
  expect_length(pars, 14)
})

test_that("parameters are correct", {
  pars <- DAISIE_create_sim_pars()
  expect_true(are_DAISIE_create_sim_pars(pars))
})

test_that("DAISIE_create_sim_pars throws error and prints correct message
          when pars[4] == 0 && nonoceanic_pars[1] == 0", {
  expect_error(
  DAISIE_create_sim_pars(pars = c(1, 1, 20, 0, 0.1),
                          nonoceanic_pars = c(0,0)),
  "Immigration rate is zero with no initial species.")
          })
