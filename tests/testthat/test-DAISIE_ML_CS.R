context("DAISIE_ML_CS")

test_that("multi-rate DAISIE_ML_CS produces correct output", {
  skip("test takes too long atm")
  utils::data(Galapagos_datalist)
  relaxed_clado <- DAISIE_ML_CS(datalist = Galapagos_datalist,
                                initparsopt = c(2, 2.7, 20, 0.009, 1.01),
                                idparsopt = 1:5,
                                parsfix = NULL,
                                idparsfix = NULL,
                                ddmodel = 11,
                                CS_version = 2.2)
})

test_that("multi-rate DAISIE_ML_CS converges to constant rate", {
  skip("WIP")
  utils::data(Galapagos_datalist)
  CR <- DAISIE_ML_CS(datalist = Galapagos_datalist,
                     initparsopt = c(2.5, 2.7, 20, 0.009, 1.01),
                     idparsopt = 1:5,
                     parsfix = NULL,
                     idparsfix = NULL,
                     ddmodel = 11,
                     CS_version = 1)
  relaxed_clado <- DAISIE_ML_CS(datalist = Galapagos_datalist,
                                initparsopt = c(2, 2.7, 20, 0.009, 1.01),
                                idparsopt = 1:5,
                                parsfix = NULL,
                                idparsfix = NULL,
                                ddmodel = 11,
                                CS_version = -2.2)
})
