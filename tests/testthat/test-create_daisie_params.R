test_that("use", {
  expect_silent(
    create_daisie_params(
       time = 4,
       M = 1,
       pars = c(2.5, 2.6, Inf, 0.01, 1.0),
       replicates = 1
    )
  )
})
