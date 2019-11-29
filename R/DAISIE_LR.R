#' Bootstrap likelihood ratio test
#'
#' @param datalist stub
#' @param datatype stub
#' @param initparsopt_dd stub
#' @param initparsopt_di stub
#' @param idparsnoshift stub
#' @param idparsmat stub
#' @param endmc stub
#' @param seed stub
#' @param alpha stub
#' @param res stub
#' @param ddmodel stub
#' @param cond stub
#' @param island_ontogeny stub
#' @param eqmodel stub
#' @param x_E stub
#' @param x_I stub
#' @param tol stub
#' @param maxiter stub
#' @param methode stub
#' @param optimmethod stub
#' @param CS_version stub
#' @param verbose stub
#' @param tolint stub
#'
#' @return stub
#' @export
#'
#' @examples stub
DAISIE_LR <- function(
  datalist,
  datatype = "single",
  initparsopt_dd,
  initparsopt_di,
  idparsnoshift = 6:10,
  idparsmat = NULL,
  endmc = 1000,
  seed = 17,
  alpha = 0.05,
  res = 100,
  ddmodel = 0,
  cond = 0,
  island_ontogeny = NA,
  eqmodel = 0,
  x_E = 0.95,
  x_I = 0.98,
  tol = c(1e-04, 1e-05, 1e-07),
  maxiter = 2000,
  methode = "lsodes",
  optimmethod = "subplex",
  CS_version = 1,
  verbose = TRUE,
  tolint = c(1E-16, 1E-10)) {
  if (verbose == TRUE) {
    cat("Estimating parameters under the diversity-independent model \n")
  }
  #estimate under DI (step 1)
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
  #estimate under DD (step 1)
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
  #calculate likelihood ratio (step 1)
  likelihood_ratio_zero <- init_dd_ml$loglik - init_di_ml$loglik

  init_di_pars <- as.numeric(c(init_di_ml[1:2], init_di_ml[4:5]))
  init_dd_pars <- as.numeric(c(init_dd_ml[1:5]))
  di_sims <- list()
  dd_sims <- list()
  time <- datalist[[1]]$island_age
  M <- datalist[[1]]$not_present + (length(datalist) - 1)
  #simulate under DI (step 2)
  cat("\nSimulating under DI \n")
  for (mc in 1:endmc) {
    di_sims[[mc]] <- DAISIE_sim(time = time,
                          M = M,
                          pars = init_di_pars,
                          replicates = 1,
                          plot_sims = FALSE,
                          verbose = FALSE)
  }
   #simulate under DD (step 6)
   cat("\nSimulating under DD \n")
    for (me in 1:endmc) {
      dd_sims[[mc]] <- DAISIE_sim(time = time,
                          M = M,
                          pars = init_dd_pars,
                          replicates = 1,
                          plot_sims = FALSE,
                          verbose = FALSE)
    }
  di_ml_est_di <- list()
  dd_ml_est_di <- list()
   #esimate under DI and DD from DI simulations (step 3)
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
  }
  di_ml_est_dd <- list()
  dd_ml_est_dd <- list()
  #estimate under DI and DD from DD simulations (step 7)
  cat('\nPerforming bootstrap to determine power ...\n')
  for(mc in 1:endmc) {
    cat('\nAnalyzing simulation:',mc,'\n')
    di_ml_est_dd[[mc]] <- DAISIE_ML_CS(datalist = dd_sims[[mc]],
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

    dd_ml_est_dd[[mc]] <- DAISIE_ML_CS(datalist = dd_sims[[mc]],
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




    LR = pmax(0,maxLLDD - outCR$loglik)



}
}

