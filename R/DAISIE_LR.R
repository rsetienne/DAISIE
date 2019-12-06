DAISIE_LR <- function(datalist,
                      datatype = "single",
                      initparsopt_di,
                      initparsopt_dd,
                      idparsopt,
                      parsfix,
                      idparsfix,
                      idparsnoshift = 6:10,
                      idparsmat = NULL,
                      res = 100,
                      ddmodel = 0,
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
                      verbose = TRUE,
                      tolint = c(1E-16, 1E-10),
                      endmc
                      ) {
  #Step 1
  if (verbose == TRUE) {
    cat("Estimating parameters under the diversity-independent model \n")
  }

  init_di_ml <- DAISIE_ML_CS(datalist = datalist,
                             datatype = datatype,
                             initparsopt = initparsopt_di,
                             idparsopt = c(1, 2, 4, 5),
                             parsfix = Inf,
                             idparsfix = 3,
                             idparsnoshift = idparsnoshift,
                             idparsmat = idparsmat,
                             res = res,
                             ddmodel = ddmodel, #should this be set to zero or does it not matter when K is set to inf?
                             cond = cond,
                             island_ontogeny = island_ontogeny,
                             eqmodel = eqmodel,
                             x_E = x_E,
                             x_I = x_I,
                             tol = tol,
                             maxiter = maxiter,
                             methode = methode,
                             optimmethod = optimmethod,
                             CS_version = CS_version,
                             verbose = verbose,
                             tolint = tolint)
  if (verbose == TRUE) {
    cat("Estimating parameters under the diversity-dependent model \n")
  }

  init_dd_ml <- DAISIE_ML_CS(datalist = datalist,
                             datatype = datatype,
                             initparsopt = initparsopt_dd,
                             idparsopt = 1:5,
                             parsfix = NULL,
                             idparsfix = NULL,
                             idparsnoshift = idparsnoshift,
                             idparsmat = idparsmat,
                             res = res,
                             ddmodel = ddmodel,
                             cond = cond,
                             island_ontogeny = island_ontogeny,
                             eqmodel = eqmodel,
                             x_E = x_E,
                             x_I = x_I,tol = tol,
                             maxiter = maxiter,
                             methode = methode,
                             optimmethod = optimmethod,
                             CS_version = CS_version,
                             verbose = verbose,
                             tolint = tolint)

  likelihood_ratio_zero <- init_dd_ml$loglik - init_di_ml$loglik

  #Step 2
  di_sims <- list()
  time <- datalist[[1]]$island_age
  M <- datalist[[1]]$not_present + (length(datalist) - 1)
  init_di_pars <- as.numeric(c(init_di_ml[1:2], Inf, init_di_ml[4:5]))
  init_dd_pars <- as.numeric(c(init_dd_ml[1:5]))
print(init_di_pars)
print(init_dd_pars)
  cat("\nSimulating under DI \n")
  for (mc in 1:endmc) {
    di_sims[[mc]] <- DAISIE_sim(time = time,
                                M = M,
                                pars = init_di_pars,
                                replicates = 1,
                                plot_sims = FALSE,
                                verbose = FALSE)
  }
  #Step 3
  di_ml_est_di <- list()
  dd_ml_est_di <- list()
  likelihood_ratio_di <- c()
print(di_sims[[1]])
  cat('\nPerforming bootstrap to determine critical LR ...\n')
  for(mc in 1:endmc) {
    cat('\nAnalyzing simulation:',mc,'\n')
    di_ml_est_di[[mc]] <- DAISIE_ML_CS(datalist = di_sims[[mc]],
                                       datatype = datatype,
                                       initparsopt = init_di_pars,
                                       idparsopt = c(1, 2, 4, 5),
                                       parsfix = Inf,
                                       idparsfix = 3,
                                       idparsnoshift = idparsnoshift,
                                       idparsmat = idparsmat,
                                       res = res,
                                       ddmodel = ddmodel, #should this be set to zero or does it not matter when K is set to inf?
                                       cond = cond,
                                       island_ontogeny = island_ontogeny,
                                       eqmodel = eqmodel,
                                       x_E = x_E,
                                       x_I = x_I,
                                       tol = tol,
                                       maxiter = maxiter,
                                       methode = methode,
                                       optimmethod = optimmethod,
                                       CS_version = CS_version,
                                       verbose = verbose,
                                       tolint = tolint)

    dd_ml_est_di[[mc]] <- DAISIE_ML_CS(datalist = di_sims[[mc]],
                                       datatype = datatype,
                                       initparsopt = init_dd_pars,
                                       idparsopt = 1:5,
                                       parsfix = NULL,
                                       idparsfix = NULL,
                                       idparsnoshift = idparsnoshift,
                                       idparsmat = idparsmat,
                                       res = res,
                                       ddmodel = ddmodel,
                                       cond = cond,
                                       island_ontogeny = island_ontogeny,
                                       eqmodel = eqmodel,
                                       x_E = x_E,
                                       x_I = x_I,tol = tol,
                                       maxiter = maxiter,
                                       methode = methode,
                                       optimmethod = optimmethod,
                                       CS_version = CS_version,
                                       verbose = verbose,
                                       tolint = tolint)

    likelihood_ratio_di[mc] <- dd_ml_est_di$loglik - di_ml_est_di$loglik
  }
  return(list(di_ml_est_di = di_ml_est_di,
              dd_ml_est_di = dd_ml_est_di,
              likelihood_ratio_di = likelihood_ratio_di))
}
