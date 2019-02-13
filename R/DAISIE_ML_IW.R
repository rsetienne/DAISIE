#' Title
#'
#' @inherit DAISIE_loglik_all_choosepar
#'
#' @export
DAISIE_loglik_IW_choosepar = function(
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
  flush.console()
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
  flush.console()
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

