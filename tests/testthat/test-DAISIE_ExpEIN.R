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

test_that("DAISIE_ExpEIN and DAISIE_ExpEIN2 give the same answer for K = Inf", {
  tvec <- 10
  pars <-  c(0.3,0.1,Inf,0.006,0)
  M <- 1000
  exp1 <- DAISIE_ExpEIN(tvec = tvec,pars = pars, M = M)
  exp2 <- DAISIE_ExpEIN2(tvec = tvec,pars = pars, M = M, res = 2000)
  testthat::expect_equal(exp1,exp2)
})

test_that("DAISIE_ExpEIN2 and DAISIE_margprobdist2 give the same expectation", {
  tvec <- c(0.1,10)
  pars <-  c(0.3,0.1,20,0.006,0.1)
  M <- 1000
  exp2 <- DAISIE_ExpEIN2(tvec = tvec,pars = pars, M = M)
  mp <- DAISIE_margprobdist2(tvec = tvec,pars = pars, M = M)
  nil2resmin1 <- matrix(rep(0:(ncol(mp$probsE) - 1), each = 2), nrow = 2, byrow = FALSE)
  nil2M2 <- matrix(rep(0:(ncol(mp$probsI) - 1), each = 2), nrow = 2, byrow = FALSE)
  End <- rowSums(mp$probsE * nil2resmin1)
  Imm <- rowSums(mp$probsI * nil2M2)
  All <- rowSums(mp$probsN * nil2resmin1)
  testthat::expect_equal(exp2$ExpE,End)
  testthat::expect_equal(exp2$ExpI,Imm)
  testthat::expect_equal(exp2$ExpN,All)

  tvec <- c(0.1,10)
  pars <-  c(0.3,0.1,20,0.006,0.1)
  M <- 1000
  initEI <- rbind(c(1,0),c(2,0),c(0,1))
  exp2 <- DAISIE_ExpEIN2(tvec = tvec,pars = pars, M = M, initEI = initEI)
  mp <- DAISIE_margprobdist2(tvec = tvec,pars = pars, M = M, initEI = initEI)
  nil2resmin1 <- matrix(rep(0:(ncol(mp$probsE) - 1), each = 2), nrow = 2, byrow = FALSE)
  nil2M2 <- matrix(rep(0:(ncol(mp$probsI) - 1), each = 2), nrow = 2, byrow = FALSE)
  End <- rowSums(mp$probsE * nil2resmin1)
  Imm <- rowSums(mp$probsI * nil2M2)
  All <- rowSums(mp$probsN * nil2resmin1)
  testthat::expect_equal(exp2$ExpE,End)
  testthat::expect_equal(exp2$ExpI,Imm)
  testthat::expect_equal(exp2$ExpN,All)
})
