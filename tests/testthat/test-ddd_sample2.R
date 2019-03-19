context("ddd_sample2")

test_that("DDD::sample with zero probability gives different series", {
  
  n <- 1000
  set.seed(42)
  event_1 <- DDD::sample2(1:4, size = n, replace = TRUE, prob = c(1, 1, 1, 0))
  testit::assert(sum(event_1 == 4) == 0)
  set.seed(42)
  event_2 <- DDD::sample2(1:3, size = n, replace = TRUE, prob = c(1, 1, 1))
  expect_false(all(event_1 == event_2))
  # Tip: use DAISIE::rng_respecting_sample  
})

test_that("base sample with zero probablity gives different series", {
  
  n <- 1000
  set.seed(42)
  event_1 <- sample(1:4, size = n, replace = TRUE, prob = c(1.0, 1.0, 1.0, 0.0))
  testit::assert(sum(event_1 == 4) == 0)
  set.seed(42)
  event_2 <- sample(1:3, size = n, replace = TRUE, prob = c(1.0, 1.0, 1.0))
  expect_false(all(event_1 == event_2))
  # Tip: use DAISIE::rng_respecting_sample  
})
