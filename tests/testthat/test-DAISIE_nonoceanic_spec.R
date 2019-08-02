context("nonoceanic_spec")

test_that("output is silent", {
  expect_silent(DAISIE_nonoceanic_spec(prob_samp = 0.5,
                                prob_nonend = 0.5,
                                mainland_n = 100))
})

test_that("output is a list of three vectors", {
  native_spec <- DAISIE_nonoceanic_spec(prob_samp = 0.1,
                                 prob_nonend = 0.9,
                                 mainland_n = 1000)
  expect_true(class(native_spec) == "list")
})

test_that("native species sampled when probability of sampling is non-zero", {
  native_spec <- DAISIE_nonoceanic_spec(prob_samp = 0.1,
                                 prob_nonend = 0.9,
                                 mainland_n = 1000)
  expect_true(is.list(native_spec))
  expect_true(is.vector(native_spec[[1]]))
  expect_true(is.numeric(native_spec[[1]]))
  expect_true(is.vector(native_spec[[2]]))
  expect_true(is.numeric(native_spec[[2]]))
  expect_true(is.vector(native_spec[[3]]))
  expect_true(is.numeric(native_spec[[3]]))
  expect_gt(length(native_spec[[1]]), 0)
  expect_gt(length(native_spec[[2]]), 0)
})

test_that("no native species are sampled with zero probability of sampling", {
  prob_samp <- 0.0
  prob_nonend <- 0.9
  mainland_n <- 1000
  native_spec <- DAISIE_nonoceanic_spec(prob_samp = prob_samp,
                                 prob_nonend = prob_nonend,
                                 mainland_n = mainland_n)
  expect_true(length(native_spec[[1]]) == 0)
  expect_true(length(native_spec[[2]]) == 0)
  expect_equal(length(native_spec[[3]]), mainland_n)
})

test_that("correct number of species are sampled with seed", {
  set.seed(17)
  spec <- DAISIE_nonoceanic_spec(prob_samp = 0.5,
                          prob_nonend = 0.5,
                          mainland_n = 100)
  expect_equivalent(spec[[1]], c(51, 33, 34,  8, 45,  4,
                                 74, 85, 31, 75, 49, 21,
                                 55, 92, 39, 81, 61, 41,
                                 58, 24, 13, 26, 72, 42))
  expect_equivalent(spec[[2]], c(80,  5, 44,  3, 12, 25,
                                 40, 17, 84,  1, 22, 79,
                                 99, 16,  9, 78, 83, 14,
                                 50, 18, 64, 20, 70, 69,
                                 53, 28, 67, 93, 73,  7,
                                 95, 30))
  expect_equivalent(spec[[3]], c(2, 4, 6, 8, 10, 11, 13,
                                 15, 19, 21, 23, 24, 26,
                                 27, 29, 31, 32, 33, 34,
                                 35, 36, 37, 38, 39, 41,
                                 42, 43, 45, 46, 47, 48,
                                 49, 51, 52, 54, 55, 56,
                                 57, 58, 59, 60, 61, 62,
                                 63, 65, 66, 68, 71, 72,
                                 74, 75, 76, 77, 81, 82,
                                 85, 86, 87, 88, 89, 90,
                                 91, 92, 94, 96, 97, 98,
                                 100))
})