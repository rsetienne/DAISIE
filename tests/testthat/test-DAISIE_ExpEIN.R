test_that("use", {
  testthat::expect_silent(DAISIE_ExpEIN(
    tvec = 4,
    pars = c(0.5, 0.1, Inf, 0.01, 0.4),
    M = 1000
  ))
})

test_that("output is named list of length 3", {
  ExpEIN_out <- DAISIE_ExpEIN(
    tvec = 4,
    pars = c(0.5, 0.1, Inf, 0.01, 0.4),
    M = 1000
  )
  testthat::expect_true(
    is.list(ExpEIN_out)
  )
  testthat::expect_length(
    ExpEIN_out, 3
  )
  testthat::expect_equal(
    names(ExpEIN_out), c("ExpE", "ExpI", "ExpN")
  )
})

test_that("use with type2", {
  ExpEIN_out <- DAISIE_ExpEIN(
    tvec = 4,
    pars = c(0.5, 0.1, Inf, 0.01, 0.4, 0.7, 0.2, Inf, 0.05, 0.1, 0.1),
    M = 1000
  )
  testthat::expect_true(
    is.list(ExpEIN_out)
  )
  testthat::expect_length(
    ExpEIN_out, 3
  )
  testthat::expect_equal(
    names(ExpEIN_out), c("ExpE", "ExpI", "ExpN")
  )
})

test_that("use with tvec == Inf", {
  ExpEIN_out <- DAISIE_ExpEIN(
    tvec = Inf,
    pars = c(0.5, 0.1, Inf, 0.01, 0.4),
    M = 1000
  )
  testthat::expect_true(
    is.list(ExpEIN_out)
  )
  testthat::expect_length(
    ExpEIN_out, 3
  )
  testthat::expect_equal(
    names(ExpEIN_out), c("ExpE", "ExpI", "ExpN")
  )
})

test_that("abuse", {
  testthat::expect_error(DAISIE_ExpEIN(
    tvec = 4,
    pars = "nonsense",
    M = 1000
  ))
})
