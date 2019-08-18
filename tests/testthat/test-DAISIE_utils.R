context("is.odd")

test_that("is.odd correctly identifies odd numbers", {
  x <- 5
  expect_silent(DAISIE:::is.odd(x))
  expect_true(DAISIE:::is.odd(x))
})

test_that("is.odd correctly identifies even numbers", {
  x <- 6
  expect_silent(DAISIE:::is.odd(x))
  expect_false((DAISIE:::is.odd(x)))
})

test_that("is.odd gives an error for a non-numeric vector", {
  x <- "string"
  expect_error(DAISIE:::is.odd(x), "'x' must be a single numeric")
})


test_that("is.odd gives an error for a vector of length > 1", {
  x <- c(1, 1, 1)
  expect_error(DAISIE:::is.odd(x))
})

context("translate_island_ontogeny")

test_that("constant ontogeny outputs a 0", {
  expect_silent(translate_island_ontogeny(0))
  expect_silent(translate_island_ontogeny("const"))
  expect_equal(translate_island_ontogeny(0), 0)
  expect_equal(translate_island_ontogeny("const"), 0)
})

test_that("linear ontogeny outputs a 1", {
  expect_silent(translate_island_ontogeny(1))
  expect_silent(translate_island_ontogeny("linear"))
  expect_equal(translate_island_ontogeny(1), 1)
  expect_equal(translate_island_ontogeny("linear"), 1)
})

test_that("beta ontogeny outputs a 2", {
  expect_silent(translate_island_ontogeny(2))
  expect_silent(translate_island_ontogeny("beta"))
  expect_equal(translate_island_ontogeny(2), 2)
  expect_equal(translate_island_ontogeny("beta"), 2)
})