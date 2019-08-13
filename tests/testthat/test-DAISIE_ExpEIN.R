context("test-DAISIE_ExpEIN")

test_that("use", {
  expect_silent(DAISIE_ExpEIN(
    t = 4,
    pars = c(0.5,0.1,Inf,0.01,0.4),
    M = 1000
  ))
})

test_that("output is named list of length 3", {
  
  ExpEIN_out <- DAISIE_ExpEIN(
    t = 4,
    pars = c(0.5,0.1,Inf,0.01,0.4),
    M = 1000
  )
  expect_true(
    is.list(ExpEIN_out)
  )
  expect_length(
    ExpEIN_out, 3
  )
  
  expect_equal(
    names(ExpEIN_out), c("ExpE", "ExpI", "ExpN")    
  )
})

test_that("abuse", {
  expect_error(DAISIE_ExpEIN(
    t = 4,
    pars = "nonsense",
    M = 1000
  ))
})

