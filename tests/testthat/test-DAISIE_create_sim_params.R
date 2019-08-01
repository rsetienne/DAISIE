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