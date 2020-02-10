#' @name DAISIE_ML
#' @aliases DAISIE_ML_CS DAISIE_ML
#' @title Maximization of the loglikelihood under the DAISIE model with clade-specific
#' diversity-dependence
#' @description This function computes the maximum likelihood estimates of the parameters of
#' the DAISIE model with clade-specific diversity-dependence for data from
#' lineages colonizing an island. It also outputs the corresponding
#' loglikelihood that can be used in model comparisons.
#' The result of sort(c(idparsopt, idparsfix, idparsnoshift)) should be
#' identical to c(1:10). If not, an error is reported that the input is
#' incoherent. The same happens when the length of initparsopt is different
#' from the length of idparsopt, and the length of parsfix is different from
#' the length of idparsfix.\cr Including the 11th parameter (p_f) in either
#' idparsopt or idparsfix (and therefore initparsopt or parsfix) is optional.
#' If this parameter is not specified, then the information in the data is
#' used, otherwise the information in the data is overruled.
#'
#' @inheritParams default_params_doc
#'
#' @return The output is a dataframe containing estimated parameters and
#' maximum loglikelihood.  \item{lambda_c}{ gives the maximum likelihood
#' estimate of lambda^c, the rate of cladogenesis} \item{mu}{ gives the maximum
#' likelihood estimate of mu, the extinction rate} \item{K}{ gives the maximum
#' likelihood estimate of K, the carrying-capacity} \item{gamma}{ gives the
#' maximum likelihood estimate of gamma, the immigration rate }
#' \item{lambda_a}{ gives the maximum likelihood estimate of lambda^a, the rate
#' of anagenesis} \item{lambda_c2}{ gives the maximum likelihood estimate of
#' lambda^c2, the rate of cladogenesis for the optional second group of
#' species} \item{mu2}{ gives the maximum likelihood estimate of mu2, the
#' extinction rate for the optional second group of species} \item{K2}{ gives
#' the maximum likelihood estimate of K2, the carrying-capacity for the
#' optional second group of species} \item{gamma2}{ gives the maximum
#' likelihood estimate of gamma2, the immigration rate for the optional second
#' group of species} \item{lambda_a2}{ gives the maximum likelihood estimate of
#' lambda^a2, the rate of anagenesis for the optional second group of species}
#' \item{loglik}{ gives the maximum loglikelihood} \item{df}{ gives the number
#' of estimated parameters, i.e. degrees of feedom} \item{conv}{ gives a
#' message on convergence of optimization; conv = 0 means convergence}
#' @author Rampal S. Etienne
#' @seealso \code{\link{DAISIE_loglik_all}},
#' \code{\link{DAISIE_sim_constant_rate}},
#' \code{\link{DAISIE_sim_time_dependent}},
#' \code{\link{DAISIE_sim_constant_rate_shift}}
#' @references Valente, L.M., A.B. Phillimore and R.S. Etienne (2015).
#' Equilibrium and non-equilibrium dynamics simultaneously operate in the
#' Galapagos islands. Ecology Letters 18: 844-852. <DOI:10.1111/ele.12461>.
#' @keywords models
#' @examples
#'
#' cat("
#' ### When all species have the same rates, and we want to optimize all 5 parameters,
#' # we use:
#'
#' utils::data(Galapagos_datalist)
#' DAISIE_ML(
#'    datalist = Galapagos_datalist,
#'    initparsopt = c(2.5,2.7,20,0.009,1.01),
#'    ddmodel = 11,
#'    idparsopt = 1:5,
#'    parsfix = NULL,
#'    idparsfix = NULL
#' )
#'
#' ### When all species have the same rates, and we want to optimize all parameters
#' # except K (which we set equal to Inf), we use:
#'
#' utils::data(Galapagos_datalist)
#' DAISIE_ML(
#'    datalist = Galapagos_datalist,
#'    initparsopt = c(2.5,2.7,0.009,1.01),
#'    idparsopt = c(1,2,4,5),
#'    parsfix = Inf,
#'    idparsfix = 3
#'    )
#'
#' ### When all species have the same rates except that the finches have a different
#' # rate of cladogenesis, and we want to optimize all parameters except K (which we
#' # set equal to Inf), fixing the proportion of finch-type species at 0.163, we use:
#'
#' utils::data(Galapagos_datalist_2types)
#' DAISIE_ML(
#'    datalist = Galapagos_datalist_2types,
#'    initparsopt = c(0.38,0.55,0.004,1.1,2.28),
#'    idparsopt = c(1,2,4,5,6),
#'    parsfix = c(Inf,Inf,0.163),
#'    idparsfix = c(3,8,11),
#'    idparsnoshift = c(7,9,10)
#'    )
#'
#' ### When all species have the same rates except that the finches have a different
#' # rate of cladogenesis, extinction and a different K, and we want to optimize all
#' # parameters, fixing the proportion of finch-type species at 0.163, we use:
#'
#' utils::data(Galapagos_datalist_2types)
#' DAISIE_ML(
#'    datalist = Galapagos_datalist_2types,
#'    ddmodel = 11,
#'    initparsopt = c(0.19,0.09,0.002,0.87,20,8.9,15),
#'    idparsopt = c(1,2,4,5,6,7,8),
#'    parsfix = c(Inf,0.163),
#'    idparsfix = c(3,11),
#'    idparsnoshift = c(9,10)
#'    )
#'
#'
#' ### When all species have the same rates except that the finches have a different
#' # rate of extinction, and we want to optimize all parameters except K (which we
#' # set equal to Inf), and we also# want to estimate the fraction of finch species
#' # in the mainland pool. we use:
#'
#' utils::data(Galapagos_datalist_2types)
#' DAISIE_ML(
#'    datalist = Galapagos_datalist_2types,
#'    initparsopt = c(2.48,2.7,0.009,1.01,2.25,0.163),
#'    idparsopt = c(1,2,4,5,7,11),
#'    parsfix = c(Inf,Inf),
#'    idparsfix = c(3,8),
#'    idparsnoshift = c(6,9,10)
#'    )
#'
#' ### When we have two islands with the same rates except for immigration and anagenesis rate,
#' # and we want to optimize all parameters, we use:
#'
#' utils::data(Galapagos_datalist)
#' DAISIE_ML(
#'    datalist = list(Galapagos_datalist,Galapagos_datalist),
#'    datatype = 'multiple',
#'    initparsopt = c(2.5,2.7,20,0.009,1.01,0.009,1.01),
#'    idparsmat = rbind(1:5,c(1:3,6,7)),
#'    idparsopt = 1:7,
#'    parsfix = NULL,
#'    idparsfix = NULL
#' )
#'
#' ### When we consider the four Macaronesia archipelagoes and set all parameters the same
#' # except for rates of cladogenesis, extinction and immigration for Canary Islands,
#' # rate of cladogenesis is fixed to 0 for the other archipelagoes,
#' # diversity-dependence is assumed to be absent
#' # and we want to optimize all parameters, we use:
#'
#' utils::data(Macaronesia_datalist)
#' DAISIE_ML(
#'    datalist = Macaronesia_datalist,
#'    datatype = 'multiple',
#'    initparsopt = c(1.053151832,0.052148979,0.512939011,0.133766934,0.152763179),
#'    idparsmat = rbind(1:5,c(6,2,3,7,5),1:5,1:5),
#'    idparsopt = c(2,4,5,6,7),
#'    parsfix = c(0,Inf),
#'    idparsfix = c(1,3)
#' )
#'
#' ")
#'
#' @export DAISIE_ML_CS
#' @export DAISIE_ML
DAISIE_ML_CS <- DAISIE_ML <- function(
     datalist,
     datatype = "single",
     initparsopt,
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
     verbose = 0,
     tolint = c(1E-16, 1E-10)
     ) {
  if (datatype == "single") {
     if (is.na(island_ontogeny)) {
       out <- DAISIE_ML1(datalist = datalist,
                        initparsopt = initparsopt,
                        idparsopt = idparsopt,
                        parsfix = parsfix,
                        idparsfix = idparsfix,
                        idparsnoshift = idparsnoshift,
                        res = res,
                        ddmodel = ddmodel,
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
                        verbose = 0,
                        tolint = tolint)
     } else {
       out <- DAISIE_ML3(datalist = datalist,
                        initparsopt = initparsopt,
                        idparsopt = idparsopt,
                        parsfix = parsfix,
                        idparsfix = idparsfix,
                        res = res,
                        ddmodel = ddmodel,
                        cond = cond,
                        island_ontogeny = island_ontogeny,
                        tol = tol,
                        maxiter = maxiter,
                        methode = methode,
                        optimmethod = optimmethod,
                        CS_version = CS_version,
                        verbose = 0,
                        tolint = tolint)
     }
  } else {
     out <- DAISIE_ML2(datalist = datalist,
                      initparsopt = initparsopt,
                      idparsopt = idparsopt,
                      parsfix = parsfix,
                      idparsfix = idparsfix,
                      idparsmat = idparsmat,
                      res = res,
                      ddmodel = ddmodel,
                      cond = cond,
                      tol = tol,
                      maxiter = maxiter,
                      methode = methode,
                      optimmethod = optimmethod,
                      verbose = 0,
                      tolint = tolint)
  }
  return(out)
}
