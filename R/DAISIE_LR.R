DAISIE_LR <- function(datalist,
                      initparsopt,
                      dmodel = 11,
                      idparsopt = 1:5,
                      parsfix = NULL,
                      idparsfix = NULL,
                      seed = 17,
                      endmc = 1000,
                      alpha = 0.05,
                      datatype = "single",
                      idparsnoshift = 6:10,
                      idparsmat = NULL,
                      res = 100,
                      cond = 0,
                      island_ontogeny = NA,
                      eqmodel = 0,
                      x_E = 0.95,
                      x_I = 0.98,
                      tol = c(1e-04, 1e-05, 1e-07),
                      maxiter = 1000 * round((1.25) ^ length(idparsopt)),
                      methode = "lsodes",
                      optimmethod = "subplex",
                      CS_version = 1,
                      verbose = 0,
                      tolint = c(1E-16, 1E-10)) {

  #DI ML model
  init_di_ml <- DAISIE_ML_CS(datalist =,
                       datatype = ,
                       initparsopt = ,
                       idparsopt = ,
                       parsfix = ,
                       idparsfix = ,
                       idparsnoshift = ,
                       idparsmat = ,
                       res = ,
                       ddmodel = ,
                       cond = ,
                       island_ontogeny = ,
                       eqmodel = ,
                       x_E = ,
                       x_I = ,
                       tol = ,
                       maxiter =,
                       methode = ,
                       optimmethod = ,
                       CS_version = ,
                       verbose = ,
                       tolint = )

  print("\nEstimating parameters under the diversity-dependent model ...\n")

  #DD ML model
  init_dd_ml <- DAISIE_ML_CS(datalist =,
                       datatype = ,
                       initparsopt = ,
                       idparsopt = ,
                       parsfix = ,
                       idparsfix = ,
                       idparsnoshift = ,
                       idparsmat = ,
                       res = ,
                       ddmodel = ,
                       cond = ,
                       island_ontogeny = ,
                       eqmodel = ,
                       x_E = ,
                       x_I = ,
                       tol = ,
                       maxiter =,
                       methode = ,
                       optimmethod = ,
                       CS_version = ,
                       verbose = ,
                       tolint = )

  likelihood_ratio_zero = init_dd_ml$loglik - init_di_ml$loglik

  di_pars <- as.numeric(outCRO[1:2])
  dd_pars <- as.numeric(outDDO[1:3])
  di_tree <- list()
  dd_tree <- list()

  cat('\nSimulating trees under CR and DD models ...\n')

  for(mc in 1:endmc) {
    di_tree[[mc]] = DAISIE_sim(time = ,
                               M = ,
                               pars = di_pars,
                               replicates = )

    dd_tree[[mc]] = DAISIE_sim(time = ,
                               M = ,
                               pars = di_pars,
                               replicates = )

  }
}
