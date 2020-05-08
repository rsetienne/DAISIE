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
  cpus <- 1
  relaxed_clado_1 <- system.time(DAISIE_ML_CS(datalist = Galapagos_datalist,
                                              initparsopt = c(2, 2.7, 20, 0.009, 1.01),
                                              idparsopt = 1:5,
                                              parsfix = NULL,
                                              idparsfix = NULL,
                                              ddmodel = 11,
                                              CS_version = list(choice = 2,
                                                                pick_parameter = 'lambda^c',
                                                                distribution = 'gamma',
                                                                sd_par = 2,
                                                                number_of_cores = cpus)))

  cpus <- 4
  cl <- initiate_cluster(cpus)
  relaxed_clado_2 <- system.time(DAISIE_ML_CS(datalist = Galapagos_datalist,
                                  initparsopt = c(2, 2.7, 20, 0.009, 1.01),
                                  idparsopt = 1:5,
                                  parsfix = NULL,
                                  idparsfix = NULL,
                                  ddmodel = 11,
                                  CS_version = list(choice = 2,
                                                 pick_parameter = 'lambda^c',
                                                 distribution = 'gamma',
                                                 sd_par = 2,
                                                 number_of_cores = cpus)))
})
