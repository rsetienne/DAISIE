test_that("bootstrapping is silent", {
  skip("WIP")
  utils::data(Galapagos_datalist)
  expect_silent(DAISIE_LR(datalist = Galapagos_datalist,
                          initparsopt_dd = c(2.5,2.7,20,0.009,1.01),
                          initparsopt_di = c(2.5,2.7,0.009,1.01)))
})

test_that("calculates maximum likelihood for data for each model", {
  skip("WIP")
  expect_true()
})
