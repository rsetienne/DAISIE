context("is_numeric_list")

test_that("is_numeric_list use", {

  expect_true(DAISIE:::is_numeric_list(x = as.list(c(3))))
  expect_true(DAISIE:::is_numeric_list(x = as.list(3)))
  expect_false(DAISIE:::is_numeric_list(x = as.list(NA)))
  expect_false(DAISIE:::is_numeric_list(as.list("Not a number")))
})
