context("is.odd")

test_that("is.odd use", {
  expect_true(is.odd(1))
  expect_false(is.odd(0))
})

test_that("is.odd abuse",{
  expect_error(is.odd("nonsense"), "'x' should be a number")
  expect_error(is.odd(NA), "'x' should be a number")
  expect_error(is.odd(NULL), "'x' should be a number")
  expect_error(is.odd(1.2), "'x' should be a whole number")
})

