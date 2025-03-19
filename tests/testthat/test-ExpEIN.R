test_that("DAISIE_ExpEIN and DAISIE_ExpEIN2 give the same answer", {
  t <- 10
  pars <-  c(0.3,0.1,Inf,0.006,0)
  M <- 1000
  exp1 <- DAISIE_ExpEIN(t = t,pars = pars, M = M)
  exp2 <- DAISIE_ExpEIN2(t = t,pars = pars, M = M)
  testthat::expect_equal(exp1,exp2)
})
