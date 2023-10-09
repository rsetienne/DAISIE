test_that("Returns TRUE when correct", {
  testthat::expect_true(
    is_sea_level_input(sea_level = "const")
  )
  testthat::expect_true(
    is_sea_level_input(sea_level = "sine")
  )
})

test_that("Returns FALSE when not string", {
  testthat::expect_false(
    is_sea_level_input(sea_level = 2)
  )
})

test_that("Returns FALSE when wrong string", {
  testthat::expect_false(
    is_sea_level_input(sea_level = "nonsense")
  )
})
