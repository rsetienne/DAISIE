test_that("The parameter choice for 2type DAISIE_ML works", {
  Galapagos_datalist_2types <- NULL
  rm(Galapagos_datalist_2types)
  utils::data(Galapagos_datalist_2types, package = "DAISIE")
  set.seed(1)
  # MLE and high tolerance for speed-up
  Fit <- DAISIE_ML(
    datalist = Galapagos_datalist_2types,
    initparsopt = c(2.183336,2.517413,0.009909,1.080458,1.316296,0.001416),
    idparsopt = c(1,2,4,5,7,11),
    parsfix = c(Inf,Inf),
    idparsfix = c(3,8),
    idparsnoshift = c(6,9,10),
    res = 30, 
    tol = c(1, 1, 1),
    maxiter = 30
  )
  testthat::expect_equal(Fit$loglik, -74.7557, tol = 1E-3)
})
