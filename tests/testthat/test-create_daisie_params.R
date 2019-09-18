test_that("use", {
  time <- 3
  M <- 1
  pars <- c(2.5, 2.6, Inf, 0.01, 1.0)
  replicates <- 1
  daisie_params <- create_daisie_params(
     time = time,
     M = M,
     pars = pars,
     replicates = replicates
  )
  expect_equal(daisie_params$time, time)
  expect_equal(daisie_params$M, M)
  expect_equal(daisie_params$pars, pars)
  expect_equal(daisie_params$replicates, replicates)

})
