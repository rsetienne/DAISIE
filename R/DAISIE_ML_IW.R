DAISIE_loglik_IW_choosepar <- function(
  trparsopt,
  trparsfix,
  idparsopt,
  idparsfix,
  pars2,
  M,
  datalist,
  methode = 'ode45',
  abstolint = 1E-16,
  reltolint = 1E-14
  )
{
  trpars1 = rep(0,5)
  trpars1[idparsopt] = trparsopt
  if(length(idparsfix) != 0)
  {
    trpars1[idparsfix] = trparsfix
  }
  if(max(trpars1) > 1 | min(trpars1) < 0)
  {
    loglik = -Inf
  } else {
    pars1 = c(trpars1/(1 - trpars1),M)
    if(min(pars1) < 0)
    {
      loglik = -Inf
    } else {
      loglik = DAISIE_loglik_IW(pars1 = pars1,pars2 = pars2,datalist = datalist,methode = methode,abstolint = abstolint,reltolint = reltolint)
    }
    if(is.nan(loglik) || is.na(loglik))
    {
      cat("There are parameter values used which cause numerical problems.\n")
      loglik = -Inf
    }
  }
  return(loglik)
}



#' Maximization of the loglikelihood under the DAISIE model with island-wide
#' diversity-dependence
#' 
#' This function computes the maximum likelihood estimates of the parameters of
#' the DAISIE model with island-wide diversity-dependence for data from
#' lineages colonizing an island. It also outputs the corresponding
#' loglikelihood that can be used in model comparisons.
#' 
#' The result of sort(c(idparsopt, idparsfix)) should be identical to c(1:5).
#' If not, an error is reported that the input is incoherent. The same happens
#' when the length of initparsopt is different from the length of idparsopt,
#' and the length of parsfix is different from the length of idparsfix.\cr
#' 
#' @param datalist Data object containing information on colonisation and
#' branching times. This object can be generated using the DAISIE_dataprep
#' function, which converts a user-specified data table into a data object, but
#' the object can of course also be entered directly. It is an R list object
#' with the following elements.\cr The first element of the list has two three
#' components: \cr \cr \code{$island_age} - the island age \cr
#' \code{$not_present} - the number of mainland lineages that are not present
#' on the island \cr The remaining elements of the list each contains
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
#' @param idparsfix The ids of the parameters that should not be optimized,
#' e.g. c(1,3) if lambda^c and K should not be optimized.
#' @param parsfix The values of the parameters that should not be optimized
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
#' @param verbose Specifies whether intermediate output should be provided,
#' because optimizationmay take long. Default is 0, no output. A value of 1
#' means the parameters and loglikelihood are printed. A value of 2 means also
#' intermediate progress during loglikelihood computation is shown.
#' @param tolint Vector of two elements containing the absolute and relative
#' tolerance of the integration
#' @return The output is a dataframe containing estimated parameters and
#' maximum loglikelihood.  \item{lambda_c}{ gives the maximum likelihood
#' estimate of lambda^c, the rate of cladogenesis} \item{mu}{ gives the maximum
#' likelihood estimate of mu, the extinction rate} \item{K}{ gives the maximum
#' likelihood estimate of K, the carrying-capacity} \item{gamma}{ gives the
#' maximum likelihood estimate of gamma, the immigration rate }
#' \item{lambda_a}{ gives the maximum likelihood estimate of lambda^a, the rate
#' of anagenesis} \item{loglik}{ gives the maximum loglikelihood} \item{df}{
#' gives the number of estimated parameters, i.e. degrees of feedom}
#' \item{conv}{ gives a message on convergence of optimization; conv = 0 means
#' convergence}
#' @author Rampal S. Etienne
#' @seealso \code{\link{DAISIE_loglik_IW}}, \code{\link{DAISIE_ML_CS}}
#' \code{\link{DAISIE_sim}}
#' @references Valente, L.M., A.B. Phillimore and R.S. Etienne (2015).
#' Equilibrium and non-equilibrium dynamics simultaneously operate in the
#' Galapagos islands. Ecology Letters 18: 844-852. <DOI:10.1111/ele.12461>.
#' @keywords models
#' @export DAISIE_ML_IW
DAISIE_ML_IW = function(
  datalist,
  initparsopt,
  idparsopt,
  parsfix,
  idparsfix,
  res = 100,
  ddmodel = 11,
  cond = 0,
  tol = c(1E-4, 1E-5, 1E-7),
  maxiter = 1000 * round((1.25)^length(idparsopt)),
  methode = "ode45",
  optimmethod = 'subplex',
  verbose = 0,
  tolint = c(1E-16,1E-14)
)
{
  options(warn=-1)
  out2err = data.frame(lambda_c = NA, mu = NA,K = NA, gamma = NA, lambda_a = NA, loglik = NA, df = NA, conv = NA)
  out2err = invisible(out2err)
  if(is.null(datalist[[1]]$brts_table))
  {
    datalist = Add_brt_table(datalist)
  }
  np = datalist[[1]]$not_present
  if(is.null(np))
  {
    np = datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2
  } 
  if(is.null(np))
  {
    cat('Number of species not present is misspecified.\n')
    return(invisible(out2err))
  }
  M = length(datalist) - 1 + np

  idpars = sort(c(idparsopt,idparsfix))
  if((prod(idpars == (1:5)) != 1) || (length(initparsopt) != length(idparsopt)) || (length(parsfix) != length(idparsfix)))
  {
    cat("The parameters to be optimized and/or fixed are incoherent.\n")
    return(out2err)
  }
  if(length(idparsopt) > 11)
  {
    cat("The number of parameters to be optimized is too high.\n")
    return(out2err)
  }
  namepars = c("lambda_c","mu","K'","gamma","lambda_a")
  if(length(namepars[idparsopt]) == 0) { optstr = "nothing" } else { optstr = namepars[idparsopt] }
  cat("You are optimizing",optstr,"\n")
  if(length(namepars[idparsfix]) == 0) { fixstr = "nothing" } else { fixstr = namepars[idparsfix] }
  cat("You are fixing",fixstr,"\n")
  cat("Calculating the likelihood for the initial parameters.","\n")
  utils::flush.console()
  trparsopt = initparsopt/(1 + initparsopt)
  trparsopt[which(initparsopt == Inf)] = 1
  trparsfix = parsfix/(1 + parsfix)
  trparsfix[which(parsfix == Inf)] = 1
  pars2 = c(res,ddmodel,cond,verbose)
  optimpars = c(tol,maxiter)
  initloglik = DAISIE_loglik_IW_choosepar(trparsopt = trparsopt,trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,M = M,pars2 = pars2,datalist = datalist,methode = methode,abstolint = tolint[1],reltolint = tolint[2])
  cat("The loglikelihood for the initial parameter values is",initloglik,"\n")
  if(initloglik == -Inf)
  {
    cat("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n")
    return(out2err)
  }
  cat("Optimizing the likelihood - this may take a while.","\n")
  utils::flush.console()
  out = DDD::optimizer(optimmethod = optimmethod,optimpars = optimpars,fun = DAISIE_loglik_IW_choosepar,trparsopt = trparsopt,idparsopt = idparsopt,trparsfix = trparsfix,idparsfix = idparsfix,M = M,pars2 = pars2,datalist = datalist,methode = methode,abstolint = tolint[1],reltolint = tolint[2])
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
  MLpars1 = rep(0,5)
  MLpars1[idparsopt] = MLpars
  if(length(idparsfix) != 0) { MLpars1[idparsfix] = parsfix }
  if(MLpars1[3] > 10^7){ MLpars1[3] = Inf }
  out2 = data.frame(lambda_c = MLpars1[1], mu = MLpars1[2], K = MLpars1[3], gamma = MLpars1[4], lambda_a = MLpars1[5], loglik = ML, df = length(initparsopt), conv = unlist(out$conv))
  s1 = sprintf('Maximum likelihood parameter estimates: lambda_c: %f, mu: %f, K: %f, gamma: %f, lambda_a: %f',MLpars1[1],MLpars1[2],MLpars1[3],MLpars1[4],MLpars1[5])
  s2 = sprintf('Maximum loglikelihood: %f',ML)
  cat("\n",s1,"\n",s2,"\n")
  return(invisible(out2))
}

