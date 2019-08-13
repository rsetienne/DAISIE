context("DAISIE_tutorial")

test_that("DAISIE_tutorial opens for os == 'windows'", {
  expect_silent(DAISIE_tutorial())
})
