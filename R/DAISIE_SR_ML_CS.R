DAISIE_SR_loglik_all_choosepar = function(
  trparsopt,
  trparsfix,
  idparsopt,
  idparsfix,
  idparsnoshift,
  pars2,
  datalist,
  methode,
  CS_version,
  abstolint = 1E-16,
  reltolint = 1E-10
  )
{
   trpars1 = rep(0,11)
   trpars1[idparsopt] = trparsopt
   if(length(idparsfix) != 0)
   {
      trpars1[idparsfix] = trparsfix
   }
   if(length(idparsnoshift) != 0)
   {
      trpars1[idparsnoshift] = trpars1[idparsnoshift - 5]
   }
   if(max(trpars1) > 1 | min(trpars1) < 0)
   {
      loglik = -Inf
   } else {
      pars1 = trpars1/(1 - trpars1)
      if(min(pars1) < 0)
      {
         loglik = -Inf
      } else {
         loglik = DAISIE_SR_loglik_CS(pars1 = pars1,pars2 = pars2,datalist = datalist,methode = methode,CS_version,abstolint = abstolint,reltolint = reltolint)
      }
      if(is.nan(loglik) || is.na(loglik))
      {
         cat("There are parameter values used which cause numerical problems.\n")
         loglik = -Inf
      }
   }
   return(loglik)
}


#' @name DAISIE_SR_ML
#' @title Maximization of the loglikelihood under the DAISIE model with clade-specific
#' diversity-dependence
#' @description This function computes the maximum likelihood estimates of the parameters of
#' the DAISIE model with clade-specific diversity-dependence and a shift in
#' parameters for data from lineages colonizing an island. It also outputs the
#' corresponding loglikelihood that can be used in model comparisons.
#' 
#' The result of sort(c(idparsopt, idparsfix, idparsnoshift)) should be
#' identical to c(1:10). If not, an error is reported that the input is
#' incoherent. The same happens when the length of initparsopt is different
#' from the length of idparsopt, and the length of parsfix is different from
#' the length of idparsfix.\cr Including the 11th parameter (p_f) in either
#' idparsopt or idparsfix (and therefore initparsopt or parsfix) is optional.
#' If this parameter is not specified, then the information in the data is
#' used, otherwise the information in the data is overruled.
#' 
#' @aliases DAISIE_SR_ML_CS DAISIE_SR_ML
#' @param datalist Data object containing information on colonisation and
#' branching times. This object can be generated using the DAISIE_dataprep
#' function, which converts a user-specified data table into a data object, but
#' the object can of course also be entered directly. It is an R list object
#' with the following elements.\cr The first element of the list has two three
#' components: \cr \cr \code{$island_age} - the island age \cr Then, depending
#' on whether a distinction between types is made, we have:\cr
#' \code{$not_present} - the number of mainland lineages that are not present
#' on the island \cr \cr The remaining elements of the list each contains
#' information on a single colonist lineage on the island and has 5
#' components:\cr \cr \code{$colonist_name} - the name of the species or clade
#' that colonized the island \cr \code{$branching_times} - island age and stem
#' age of the population/species in the case of Non-endemic, Non-endemic_MaxAge
#' and Endemic anagenetic species. For cladogenetic species these should be
#' island age and branching times of the radiation including the stem age of
#' the radiation.\cr \code{$stac} - the status of the colonist \cr \cr *
#' Non_endemic_MaxAge: 1 \cr * Endemic: 2 \cr * Endemic&Non_Endemic: 3 \cr *
#' Non_endemic: 4 \cr * Endemic_MaxAge: 5 \cr \cr \code{$missing_species} -
#' number of island species that were not sampled for particular clade (only
#' applicable for endemic clades) \cr
#' @param initparsopt The initial values of the parameters that must be
#' optimized
#' @param idparsopt The ids of the parameters that must be optimized. The ids
#' are defined as follows: \cr \cr id = 1 corresponds to lambda^c (cladogenesis
#' rate) \cr id = 2 corresponds to mu (extinction rate) \cr id = 3 corresponds
#' to K (clade-level carrying capacity) \cr id = 4 corresponds to gamma
#' (immigration rate) \cr id = 5 corresponds to lambda^a (anagenesis rate) \cr
#' id = 6 corresponds to lambda^c (cladogenesis rate) after the shift \cr id =
#' 7 corresponds to mu (extinction rate) after the shift\cr id = 8 corresponds
#' to K (clade-level carrying capacity) after the shift\cr id = 9 corresponds
#' to gamma (immigration rate) after the shift\cr id = 10 corresponds to
#' lambda^a (anagenesis rate) after the shift\cr id = 11 corresponds to the
#' time of shift
#' @param idparsfix The ids of the parameters that should not be optimized,
#' e.g. c(1,3) if lambda^c and K should not be optimized.
#' @param parsfix The values of the parameters that should not be optimized
#' @param idparsnoshift The ids of the parameters that should not be different
#' before and after the shift
#' @param res Sets the maximum number of species for which a probability must
#' be computed, must be larger than the size of the largest clade
#' @param ddmodel Sets the model of diversity-dependence: \cr \cr ddmodel = 0 :
#' no diversity dependence \cr ddmodel = 1 : linear dependence in speciation
#' rate \cr ddmodel = 11: linear dependence in speciation rate and in
#' immigration rate \cr ddmodel = 2 : exponential dependence in speciation
#' rate\cr ddmodel = 21: exponential dependence in speciation rate and in
#' immigration rate\cr
#' @param cond cond = 0 : conditioning on island age \cr cond = 1 :
#' conditioning on island age and non-extinction of the island biota \cr
#' @param island_ontogeny type of island ontonogeny. If NA, then constant ontogeny is assumed
#' @param tol Sets the tolerances in the optimization. Consists of: \cr reltolx
#' = relative tolerance of parameter values in optimization \cr reltolf =
#' relative tolerance of function value in optimization \cr abstolx = absolute
#' tolerance of parameter values in optimization
#' @param maxiter Sets the maximum number of iterations in the optimization
#' @param methode Method of the ODE-solver. See package deSolve for details.
#' Default is "lsodes"
#' @param optimmethod Method used in likelihood optimization. Default is
#' "subplex" (see subplex package). Alternative is 'simplex' which was the
#' method in previous versions.
#' @param CS_version For internal testing purposes only. Default is 1, the
#' original DAISIE code.
#' @param verbose sets whether parameters and likelihood should be printed (1)
#' or not (0)
#' @param tolint Vector of two elements containing the absolute and relative
#' tolerance of the integration
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
#' @seealso \code{\link{DAISIE_loglik_all}}, \code{\link{DAISIE_sim}}
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
#' @export DAISIE_SR_ML_CS
#' @export DAISIE_SR_ML
DAISIE_SR_ML_CS <- DAISIE_SR_ML <- function(
  datalist,
  initparsopt,
  idparsopt,
  parsfix,
  idparsfix,
  idparsnoshift = 6:10,
  res = 100,
  ddmodel = 0,
  cond = 0,
  island_ontogeny = NA,
  tol = c(1E-4, 1E-5, 1E-7),
  maxiter = 1000 * round((1.25)^length(idparsopt)),
  methode = "lsodes",
  optimmethod = 'subplex',
  CS_version = 1,
  verbose = 0,
  tolint = c(1E-16,1E-10)
  )
{
# datalist = list of all data: branching times, status of clade, and numnber of missing species
# datalist[[,]][1] = list of branching times (positive, from present to past)
# - max(brts) = age of the island
# - next largest brts = stem age / time of divergence from the mainland
# The interpretation of this depends on stac (see below)
# For stac = 0, this needs to be specified only once.
# For stac = 1, this is the time since divergence from the immigrant's sister on the mainland.
# The immigrant must have immigrated at some point since then.
# For stac = 2 and stac = 3, this is the time since divergence from the mainland.
# The immigrant that established the clade on the island must have immigrated precisely at this point.
# For stac = 3, it must have reimmigrated, but only after the first immigrant had undergone speciation.
# - min(brts) = most recent branching time (only for stac = 2, or stac = 3)
# datalist[[,]][2] = list of status of the clades formed by the immigrant
#  . stac == 0 : immigrant is not present and has not formed an extant clade
# Instead of a list of zeros, here a number must be given with the number of clades having stac = 0
#  . stac == 1 : immigrant is present but has not formed an extant clade
#  . stac == 2 : immigrant is not present but has formed an extant clade
#  . stac == 3 : immigrant is present and has formed an extant clade
#  . stac == 4 : immigrant is present but has not formed an extant clade, and it is known when it immigrated.
#  . stac == 5 : immigrant is not present and has not formed an extent clade, but only an endemic species
# datalist[[,]][3] = list with number of missing species in clades for stac = 2 and stac = 3;
# for stac = 0 and stac = 1, this number equals 0.
# initparsopt, parsfix = optimized and fixed model parameters
# - pars1[1] = lac = (initial) cladogenesis rate
# - pars1[2] = mu = extinction rate
# - pars1[3] = K = maximum number of species possible in the clade
# - pars1[4] = gam = (initial) immigration rate
# - pars1[5] = laa = (initial) anagenesis rate
# - pars1[6]...pars1[10] = same as pars1[1]...pars1[5], but after the shift
# - pars1[11] = time of shift
# idparsopt, idparsfix = ids of optimized and fixed model parameters
# - res = pars2[1] = lx = length of ODE variable x
# - ddmodel = pars2[2] = diversity-dependent model,mode of diversity-dependence
#  . ddmodel == 0 : no diversity-dependence
#  . ddmodel == 1 : linear dependence in speciation rate (anagenesis and cladogenesis)
#  . ddmodel == 11 : linear dependence in speciation rate and immigration rate
#  . ddmodel == 3 : linear dependence in extinction rate
# - cond = conditioning; currently only cond = 0 is possible
#  . cond == 0 : no conditioning
#  . cond == 1 : conditioning on presence on the island

  options(warn=-1)
  out2err = data.frame(lambda_c = NA, mu = NA,K = NA, gamma = NA, lambda_a2 = NA, lambda_c2 = NA, mu2 = NA,K2 = NA, gamma2 = NA, lambda_a2 = NA, tshift = NA,loglik = NA, df = NA, conv = NA)
  out2err = invisible(out2err)
  idpars = sort(c(idparsopt,idparsfix,idparsnoshift))
  missnumspec = unlist(lapply(datalist,function(list) {list$missing_species}))
  if(CS_version != 1)
  {
    cat('This version of CS is not yet implemented\n')
    return(out2err)
  }
  if(sum(missnumspec) > (res - 1))
  {
    cat("The number of missing species is too large relative to the resolution of the ODE.\n")
    return(out2err)
  }
  if((prod(idpars == (1:11)) != 1) || (length(initparsopt) != length(idparsopt)) || (length(parsfix) != length(idparsfix)))
  {
    cat("The parameters to be optimized and/or fixed are incoherent.\n")
    return(out2err)
  }
  if(length(idparsopt) > 11)
  {
    cat("The number of parameters to be optimized is too high.\n")
    return(out2err)
  }
  namepars = c("lambda_c","mu","K","gamma","lambda_a","lambda_c2","mu2","K2","gamma2","lambda_a2","tshift")
  if(length(namepars[idparsopt]) == 0) { optstr = "nothing" } else { optstr = namepars[idparsopt] }
  cat("You are optimizing",optstr,"\n")
  if(length(namepars[idparsfix]) == 0) { fixstr = "nothing" } else { fixstr = namepars[idparsfix] }
  cat("You are fixing",fixstr,"\n")
  if(sum(idparsnoshift == (6:10)) != 5)
  {
    noshiftstring = namepars[idparsnoshift]
    cat("You are not shifting",noshiftstring,"\n")
  }
  cat("Calculating the likelihood for the initial parameters.","\n")
  utils::flush.console()
  trparsopt = initparsopt/(1 + initparsopt)
  trparsopt[which(initparsopt == Inf)] = 1
  trparsfix = parsfix/(1 + parsfix)
  trparsfix[which(parsfix == Inf)] = 1
  pars2 = c(res,ddmodel,cond,verbose,island_ontogeny,tol,maxiter)
  optimpars = c(tol,maxiter)
  initloglik = DAISIE_SR_loglik_all_choosepar(trparsopt = trparsopt,trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,idparsnoshift = idparsnoshift,pars2 = pars2,datalist = datalist,methode = methode,CS_version = CS_version,abstolint = tolint[1],reltolint = tolint[2])
  cat("The loglikelihood for the initial parameter values is",initloglik,"\n")
  if(initloglik == -Inf)
  {
    cat("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n")
    return(out2err)
  }
  cat("Optimizing the likelihood - this may take a while.","\n")
  utils::flush.console()
  out = DDD::optimizer(optimmethod = optimmethod,optimpars = optimpars,fun = DAISIE_SR_loglik_all_choosepar,trparsopt = trparsopt,idparsopt = idparsopt,trparsfix = trparsfix,idparsfix = idparsfix,idparsnoshift = idparsnoshift,pars2 = pars2,datalist = datalist,methode = methode,CS_version = CS_version,abstolint = tolint[1],reltolint = tolint[2])        
  if(out$conv != 0)
  {
    cat("Optimization has not converged. Try again with different initial values.\n")
    out2 = out2err
    out2$conv = out$conv
    return(out2)
  }
  MLtrpars = as.numeric(unlist(out$par))
  MLpars = MLtrpars/(1-MLtrpars)
  ML = as.numeric(unlist(out$fvalues))
  MLpars1 = rep(0,11)
  MLpars1[idparsopt] = MLpars
  if(length(idparsfix) != 0)
  {
    MLpars1[idparsfix] = parsfix
  }
  if(MLpars1[3] > 10^7)
  {
    MLpars1[3] = Inf
  }
  if(length(idparsnoshift) != 0)
  {
    MLpars1[idparsnoshift] = MLpars1[idparsnoshift - 5]
  }
  if(MLpars1[8] > 10^7)
  {
    MLpars1[8] = Inf
  }
  out2 = data.frame(lambda_c = MLpars1[1], mu = MLpars1[2], K = MLpars1[3], gamma = MLpars1[4], lambda_a = MLpars1[5], lambda_c2 = MLpars1[6], mu2 = MLpars1[7], K2 = MLpars1[8], gamma2 = MLpars1[9], lambda_a2 = MLpars1[10], tshift = MLpars1[11], loglik = ML, df = length(initparsopt), conv = unlist(out$conv))
  s1 = sprintf('Maximum likelihood parameter estimates: lambda_c: %f, mu: %f, K: %f, gamma: %f, lambda_a: %f, lambda_c2: %f, mu2: %f, K2: %f, gamma2: %f, lambda_a2: %f, time of shift: %f',MLpars1[1],MLpars1[2],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6],MLpars1[7],MLpars1[8],MLpars1[9],MLpars1[10],MLpars1[11])
  s2 = sprintf('Maximum loglikelihood: %f',ML)
  cat("\n",s1,"\n",s2,"\n")
  return(invisible(out2))
}
