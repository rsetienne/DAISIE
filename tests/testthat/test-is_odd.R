context("is.odd")

test_that("is.odd use", {
  skip("Not now, Issue 71, Issue #71")
  expect_true(is.odd(1))
  expect_false(is.odd(0))
})

test_that("is.odd abuse",{
  skip("Not now, Issue 71, Issue #71")
  expect_error(is.odd("nonsense"), "'x' must be numeric")
  expect_error(is.odd(NA), "'x' must be numeric")
  expect_error(is.odd(NULL), "'x' must be numeric")
  expect_error(is.odd(1.2), "'x' must be an integer")
})

