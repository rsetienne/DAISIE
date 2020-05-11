context("multi-rate integration test")
test_that("multi-rate loglik works", {
  Galapagos_datalist = NULL
  rm(Galapagos_datalist)
  utils::data(Galapagos_datalist, package = "DAISIE")
  CR <- DAISIE_ML_CS(datalist = Galapagos_datalist,
                     initparsopt = c(2.5, 2.7, 20, 0.009, 1.01),
                     idparsopt = 1:5,
                     parsfix = NULL,
                     idparsfix = NULL,
                     ddmodel = 11,
                     CS_version = 1)
  CS_version <- create_CS_version(model = "multi",
                                  pick_parameter = "cladogenesis",
                                  distribution = "gamma",
                                  sd = 2,
                                  num_cores = 1)
  relaxed_clado_1 <- system.time(DAISIE_ML_CS(datalist = Galapagos_datalist,
                                              initparsopt = c(2, 2.7, 20, 0.009, 1.01),
                                              idparsopt = 1:5,
                                              parsfix = NULL,
                                              idparsfix = NULL,
                                              ddmodel = 11,
                                              CS_version = CS_version))

  cpus <- 7
  CS_version <- create_CS_version(model = "multi",
                                  pick_parameter = "cladogenesis",
                                  distribution = "gamma",
                                  sd = 2,
                                  num_cores = cpus)
  cl <- initiate_cluster(cpus)
  relaxed_clado_2 <- system.time(DAISIE_ML_CS(datalist = Galapagos_datalist,
                                  initparsopt = c(2, 2.7, 20, 0.009, 1.01),
                                  idparsopt = 1:5,
                                  parsfix = NULL,
                                  idparsfix = NULL,
                                  ddmodel = 11,
                                  CS_version = CS_version))
})
