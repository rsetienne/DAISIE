context("DAISIE_utils")

test_that("sample_relaxed_rate produces correct output for cladogenesis", {
  pars <- c(1, 1, 20, 0.1, 1)
  relaxed_par <- "cladogenesis"
  relaxed_rate_pars <- list(mean = 1, sd = 1)
  set.seed(1)
  pars <- sample_relaxed_rate(pars = pars,
                              relaxed_par = relaxed_par,
                              relaxed_rate_pars = relaxed_rate_pars)
  expect_true(is.numeric(pars))
  expect_length(pars, 5)
  expect_equal(pars,
               c(0.155141356572, 1.0000000, 20.0000000, 0.1000000, 1.0000000))
})

test_that("sampled_relaxed_rate produces correct output for extinction", {
  pars <- c(1, 1, 20, 0.1, 1)
  relaxed_par <- "extinction"
  relaxed_rate_pars <- list(mean = 1, sd = 1)
  set.seed(1)
  pars <- sample_relaxed_rate(pars = pars,
                              relaxed_par = relaxed_par,
                              relaxed_rate_pars = relaxed_rate_pars)
  expect_true(is.numeric(pars))
  expect_length(pars, 5)
  expect_equal(pars,
               c(1.0000000, 0.155141356572, 20.0000000, 0.1000000, 1.0000000))
})

test_that("sampled_relaxed_rate produces correct output for carrying capacity", {
  pars <- c(1, 1, 20, 0.1, 1)
  relaxed_par <- "carrying_capacity"
  relaxed_rate_pars <- list(mean = 1, sd = 1)
  set.seed(1)
  pars <- sample_relaxed_rate(pars = pars,
                              relaxed_par = relaxed_par,
                              relaxed_rate_pars = relaxed_rate_pars)
  expect_true(is.numeric(pars))
  expect_length(pars, 5)
  expect_equal(pars,
               c(1.0000000, 1.0000000, 0.155141356572, 0.1000000, 1.0000000))
})

test_that("sampled_relaxed_rate produces correct output for immigration", {
  pars <- c(1, 1, 20, 0.1, 1)
  relaxed_par <- "immigration"
  relaxed_rate_pars <- list(mean = 1, sd = 1)
  set.seed(1)
  pars <- sample_relaxed_rate(pars = pars,
                              relaxed_par = relaxed_par,
                              relaxed_rate_pars = relaxed_rate_pars)
  expect_true(is.numeric(pars))
  expect_length(pars, 5)
  expect_equal(pars,
               c(1.0000000, 1.0000000, 20.0000000, 0.155141356572, 1.0000000))
})

test_that("sampled_relaxed_rate produces correct output for anagenesis", {
  pars <- c(1, 1, 20, 0.1, 1)
  relaxed_par <- "anagenesis"
  relaxed_rate_pars <- list(mean = 1, sd = 1)
  set.seed(1)
  pars <- sample_relaxed_rate(pars = pars,
                              relaxed_par = relaxed_par,
                              relaxed_rate_pars = relaxed_rate_pars)
  expect_true(is.numeric(pars))
  expect_length(pars, 5)
  expect_equal(pars,
               c(1.0000000, 1.0000000, 20.0000000, 0.1000000, 0.155141356572))
})
