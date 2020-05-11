context("DAISIE_ML_CS")


test_that("multi-rate DAISIE_ML_CS produces correct output", {
  skip("test takes too long atm")
  utils::data(Galapagos_datalist)
  cpus <- parallel::detectCores()
  CS_version <- create_CS_version(model = 2,
                                  pick_parameter = "cladogenesis",
                                  distribution = "gamma",
                                  sd = 2,
                                  num_cores = cpus)
  cl <- initiate_cluster(cpus)
  RR_clado <- system.time(DAISIE_ML_CS(datalist = Galapagos_datalist,
                                       initparsopt = c(2, 2.7, 20, 0.009, 1.01),
                                       idparsopt = 1:5,
                                       parsfix = NULL,
                                       idparsfix = NULL,
                                       ddmodel = 11,
                                       CS_version = CS_version))
  expect_true(is.numeric(RR_clado))
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

  utils::data(Galapagos_datalist)
  cpus <- 7
  CS_version <- create_CS_version(model = 2,
                                  pick_parameter = "cladogenesis",
                                  distribution = "lognormal",
                                  sd = 1,
                                  num_cores = cpus)
  RR_clado <- DAISIE_ML_CS(datalist = Galapagos_datalist,
                                initparsopt = c(2.5, 2.7, 20, 0.009, 1.01),
                                idparsopt = 1:5,
                                parsfix = NULL,
                                idparsfix = NULL,
                                ddmodel = 11,
                                CS_version = CS_version)
  expected_equal(CR, RR_clado)
})
