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
  maxiter = 1000 * round((1.25) ^ length(idparsopt)),
  methode = "lsodes",
  optimmethod = "subplex",
  CS_version = 1,
  verbose = TRUE,
  tolint = c(1E-16, 1E-10)) {
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
