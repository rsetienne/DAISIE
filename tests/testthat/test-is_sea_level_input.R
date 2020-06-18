context("is_sea_level_input")

test_that("Returns TRUE when correct", {
  expect_true(
    is_sea_level_input(sea_level = "const")
  )
  expect_true(
    is_sea_level_input(sea_level = "sine")
  )
})

test_that("Returns FALSE when not string", {
  expect_false(
    is_sea_level_input(sea_level = 2)
  )
})

test_that("Returns FALSE when wrong string", {
  expect_false(
    is_sea_level_input(sea_level = "nonsense")
  )
})
