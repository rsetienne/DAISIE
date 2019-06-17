context("rng_respecting_sample")

test_that("use, zero probability at end", {

  n <- 1000
  set.seed(42)
  draws_1 <- DAISIE:::rng_respecting_sample(
    1:4, size = n, replace = TRUE, prob = c(1.0, 1.0, 1.0, 0.0)
  )
  testit::assert(sum(draws_1 == 4) == 0)
  set.seed(42)
  draws_2 <- DAISIE:::rng_respecting_sample(
    1:3, size = n, replace = TRUE, prob = c(1.0, 1.0, 1.0)
  )
  testit::assert(sum(draws_2 == 4) == 0)
  expect_equal(draws_1, draws_2)
})

test_that("use", {

  n <- 1000
  set.seed(42)
  draws_1 <- DAISIE:::rng_respecting_sample(
    c(1, 3, 4), size = n, replace = TRUE, prob = c(1.0, 1.0, 1.0)
  )
  testit::assert(sum(draws_1 == 2) == 0)
  testit::assert(sum(draws_1 == 3) > 0)
  set.seed(42)
  draws_2 <- DAISIE:::rng_respecting_sample(
    1:4, size = n, replace = TRUE, prob = c(1.0, 0.0, 1.0, 1.0)
  )
  testit::assert(sum(draws_2 == 2) == 0)
  testit::assert(sum(draws_2 == 3) > 0)
  expect_equal(draws_1, draws_2)
})

test_that("use, two non-zero probabilities", {

  n <- 1000
  set.seed(42)
  draws_1 <- DAISIE:::rng_respecting_sample(
    c(2, 4), size = n, replace = TRUE, prob = c(1.0, 1.0)
  )
  testit::assert(sum(draws_1 == 1) == 0)
  testit::assert(sum(draws_1 == 2) > 0)
  set.seed(42)
  draws_2 <- DAISIE:::rng_respecting_sample(
    1:4, size = n, replace = TRUE, prob = c(0.0, 1.0, 0.0, 1.0)
  )
  testit::assert(sum(draws_2 == 1) == 0)
  testit::assert(sum(draws_2 == 2) > 0)
  expect_equal(draws_1, draws_2)
})
