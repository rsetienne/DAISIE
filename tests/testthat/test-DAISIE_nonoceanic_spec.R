context("DAISIE_nonoceanic_spec")

test_that("native species sampled when probability of sampling is non-zero", {
  native_spec <- DAISIE:::DAISIE_nonoceanic_spec(prob_samp = 0.1, 
                                                 prob_nonend = 0.9, 
                                                 mainland_n = 1000)
  
  expect_true(is.list(native_spec))
  expect_true(is.vector(native_spec[1]))
  expect_true(is.vector(native_spec[2]))
  expect_true(is.vector(native_spec[3]))
  expect_gt(length(native_spec[1]), 0)
  expect_gt(length(native_spec[2]), 0)
})


test_that("no native species are sampled when probability of sampling is zero", {
  native_spec <- DAISIE:::DAISIE_nonoceanic_spec(prob_samp = 0.0, 
                                                 prob_nonend = 0.9, 
                                                 mainland_n = 1000)
  
  expect_true(length(native_spec[1]) == 0)
  expect_true(length(native_spec[2]) == 0)
  expect_equal(length(native_spec[3]), mainland_n)
})