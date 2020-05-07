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
  cpus <- 4
  cpus <- min(cpus,parallel::detectCores())
  cl <- parallel::makeCluster(cpus)
  doParallel::registerDoParallel(cl, cores = cpus)
  on.exit(parallel::stopCluster(cl))
  relaxed_clado_2 <- system.time(DAISIE_ML_CS(datalist = Galapagos_datalist,
                                  initparsopt = c(2, 2.7, 20, 0.009, 1.01),
                                  idparsopt = 1:5,
                                  parsfix = NULL,
                                  idparsfix = NULL,
                                  ddmodel = 11,
                                  CS_version = 2.2))
})
