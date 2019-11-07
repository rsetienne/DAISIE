context("DAISIE_create_sim_pars")

test_that("output is silent", {
  expect_silent(DAISIE_create_sim_pars())
})

test_that("output is a list", {
  pars <- DAISIE_create_sim_pars()
  expect_true(class(pars) == "list")
  expect_length(pars, 16)
})

test_that("parameters are correct", {
  pars <- DAISIE_create_sim_pars()
  expect_true(are_DAISIE_create_sim_pars(pars))
})

test_that("DAISIE_create_sim_pars throws error and prints correct message
          when pars[4] == 0 && island_type == 'oceanic'", {
  expect_error(
  DAISIE_create_sim_pars(pars = c(1, 1, 20, 0, 0.1),
                           island_type = "oceanic"),
  "Immigration rate is zero with no initial species.")
          })

test_that("DAISIE_create_sim_pars throws warning and prints correct message
          when island_type == 'oceanic' && !is.null(nonoceanic_pars)", {
  expect_warning(
    DAISIE_create_sim_pars(island_type = "oceanic",
                             nonoceanic_pars = c(0.1, 0.9)),
    "Nonoceanic parameters have been specified with an oceanic
    island. Set nonoceanic_pars to NULL")
          })

test_that("DAISIE_create_sim_pars throws error and prints correct message
          when island_type == 'nonoceanic' && is.null(nonoceanic_pars)", {
    expect_error(
      DAISIE_create_sim_pars(island_type = "nonoceanic",
                               nonoceanic_pars = NULL),
      "Nonoceanic island has no parameters.")
          })
