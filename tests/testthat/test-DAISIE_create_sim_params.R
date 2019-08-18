context("DAISIE_create_sim_params")

test_that("output is silent", {
  expect_silent(DAISIE_create_sim_params())
})

test_that("output is a list", {
  params <- DAISIE_create_sim_params()
  expect_true(class(params) == "list")
  expect_length(params, 19)
})

test_that("parameters are correct", {
  params <- DAISIE_create_sim_params()
  expect_true(are_DAISIE_create_sim_params(params))
})

test_that("DAISIE_create_sim_params throws error and prints correct message
          when pars[4] == 0 && island_type == 'oceanic'", {
  expect_error(
  DAISIE_create_sim_params(pars = c(1,1,20,0,0.1),
                           island_type = "oceanic"),
  "Immigration rate is zero with no initial species.")
          })

test_that("DAISIE_create_sim_params throws warning and prints correct message
          when island_type == 'oceanic' && !is.null(nonoceanic)", {
  expect_warning(
    DAISIE_create_sim_params(island_type = "oceanic",
                             nonoceanic = c(0.1, 0.9)),
    "Nonoceanic parameters have been specified with an oceanic
    island. Set nonoceanic to NULL")
          })

test_that("DAISIE_create_sim_params throws error and prints correct message
          when island_type == 'nonoceanic' && is.null(nonoceanic)", {
    expect_error(
      DAISIE_create_sim_params(island_type = "nonoceanic",
                               nonoceanic = NULL),
      "Nonoceanic island has no parameters.")
          })