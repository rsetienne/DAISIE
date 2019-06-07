#' ML stub
#'
#' @inherit DAISIE_loglik_all_choosepar
#'
#' @export
DAISIE_loglik_all_choosepar3 = function(
  trparsopt,
  trparsfix,
  idparsopt,
  idparsfix,
  pars2,
  datalist,
  methode,
  CS_version = 1,
  abstolint = 1E-16,
  reltolint = 1E-10
  )
{
   trpars1 = rep(0,10)
   trpars1[idparsopt] = trparsopt
   if(length(idparsfix) != 0)
   {
      trpars1[idparsfix] = trparsfix
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
         loglik = DAISIE_loglik_all(pars1 = pars1,pars2 = pars2,datalist = datalist,methode = methode,CS_version = CS_version,abstolint = abstolint,reltolint = reltolint)
      }
      if(is.nan(loglik) || is.na(loglik))
      {
         cat("There are parameter values used which cause numerical problems.\n")
         loglik = -Inf
      }
   }
   return(loglik)
}

DAISIE_ML3 = function(
  datalist,
  initparsopt,
  idparsopt,
  parsfix,
  idparsfix,
  res = 100,
  ddmodel = 0,
  cond = 0,
  island_ontogeny,
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
# - pars1[1:4] = Apars  
# - pars1[5] = lac = (initial) cladogenesis rate
# - pars1[6:7] = extinction rate parameters
# - pars1[8] = K = maximum number of species possible in the clade
# - pars1[9] = gam = (initial) immigration rate
# - pars1[10] = laa = (initial) anagenesis rate
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
  out2err = data.frame(lambda_c = NA, mu = NA,K = NA, gamma = NA, lambda_a = NA, loglik = NA, df = NA, conv = NA)
  out2err = invisible(out2err)
  idpars = sort(c(idparsopt,idparsfix))
  missnumspec = unlist(lapply(datalist,function(list) {list$missing_species}))
  if(sum(missnumspec) > (res - 1))
  {
    cat("The number of missing species is too large relative to the resolution of the ODE.\n")
    return(out2err)
  }
  if((prod(idpars == (1:10)) != 1) || (length(initparsopt) != length(idparsopt)) || (length(parsfix) != length(idparsfix)))
  {
    cat("The parameters to be optimized and/or fixed are incoherent.\n")
    return(out2err)
  }
  if(length(idparsopt) > 10)
  {
    cat("The number of parameters to be optimized is too high.\n")
    return(out2err)
  } 
  namepars = c("Apars1","Apars2","Apars3","Apars4","lambda_c0","mu_1","mu_2","K0","gamma0","lambda_a")
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
  island_ontogeny <- translate_island_ontogeny(island_ontogeny)
  pars2 = c(res,ddmodel,cond,verbose,island_ontogeny,eqmodel = NA,tol,maxiter)
  optimpars = c(tol,maxiter)
  initloglik = DAISIE_loglik_all_choosepar3(trparsopt = trparsopt,trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,pars2 = pars2,datalist = datalist,methode = methode, CS_version = CS_version, abstolint = tolint[1], reltolint = tolint[2])
  cat("The loglikelihood for the initial parameter values is",initloglik,"\n")
  if(initloglik == -Inf)
  {
    cat("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n")
    return(out2err)
  }  
  cat("Optimizing the likelihood - this may take a while.","\n")
  utils::flush.console()
  out = DDD::optimizer(optimmethod = optimmethod,optimpars = optimpars,fun = DAISIE_loglik_all_choosepar3,trparsopt = trparsopt,idparsopt = idparsopt,trparsfix = trparsfix,idparsfix = idparsfix,pars2 = pars2,datalist = datalist,methode = methode,CS_version = CS_version,abstolint = tolint[1],reltolint = tolint[2])        
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
  MLpars1 = rep(0,10)
  MLpars1[idparsopt] = MLpars
  if(length(idparsfix) != 0) { MLpars1[idparsfix] = parsfix }
  if(MLpars1[8] > 10^7){ MLpars1[8] = Inf }
  out2 = data.frame(Apars1 = MLpars1[1], Apars2 = MLpars1[2], Apars3 = MLpars1[3], Apars4 = MLpars1[4], lambda_c0 = MLpars1[5], mu1 = MLpars1[6], mu2 = MLpars1[7], K0 = MLpars1[8], gamma0 = MLpars1[9], lambda_a = MLpars1[10],loglik = ML, df = length(initparsopt), conv = unlist(out$conv))
  s1 = sprintf('Maximum likelihood parameter estimates: Apars1: %f, Apars2: %f, Apars3: %f, Apars4: %f, lambda_c0: %f, mu1: %f, mu2: %f, K0: %f, gamma0: %f, lambda_a: %f',MLpars1[1],MLpars1[2],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6],MLpars1[7],MLpars1[8],MLpars1[9],MLpars1[10])
  s2 = sprintf('Maximum loglikelihood: %f',ML)
  cat("\n",s1,"\n",s2,"\n")
  return(invisible(out2))
}