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

test_that("abuse", {
  time <- 3
  M <- 1
  pars <- c(2.5, 2.6, Inf, 0.01, 1.0)
  replicates <- 1
  expect_silent(
    create_daisie_params(
       time = time,
       M = M,
       pars = pars,
       replicates = replicates
    )
  )
  expect_error(
    create_daisie_params(
       time = -1234567890,
       M = M,
       pars = pars,
       replicates = replicates
    ),
    "'time' must be non-zero and positive"
  )
  expect_error(
    create_daisie_params(
       time = c(1, 2, 3, 4, 5),
       M = M,
       pars = pars,
       replicates = replicates
    ),
    "'time' must be one non-zero and positive value"
  )
  expect_error(
    create_daisie_params(
       time = time,
       M = -1234567890,
       pars = pars,
       replicates = replicates
    ),
    "'M' must be non-zero and positive"
  )
  expect_error(
    create_daisie_params(
       time = time,
       M = c(1, 2, 3, 4, 5, 6),
       pars = pars,
       replicates = replicates
    ),
    "'M' must be one non-zero and positive value"
  )
  expect_error(
    create_daisie_params(
       time = time,
       M = M,
       pars = c(3, 14, 15),
       replicates = replicates
    ),
    "'pars' must have a length of at least 5"
  )
  expect_error(
    create_daisie_params(
       time = time,
       M = M,
       pars = c(-3, -1, -4, -1, -5),
       replicates = replicates
    ),
    "'pars' must be non-zero and positive"
  )
  expect_error(
    create_daisie_params(
       time = time,
       M = M,
       pars = pars,
       replicates = -12345678
    ),
    "'replicates' must be non-zero and positive"
  )

})
