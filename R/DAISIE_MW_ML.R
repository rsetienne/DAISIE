f_sigmoidal_pos <- function(d,x,d0,k)
{
  if(d0 < d)
  {
    out <- k / (1 + (d0/d)^x)
  } else
  {
    out <- k * (d/d0)^x / (1 + (d/d0)^x)
  }
  return(out)
}

f_sigmoidal_neg <- function(d,x,d0,k)
{
  out <- k - f_sigmoidal_pos(d,x,d0,k)
  return(out)
}


#' @importFrom foreach foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar%
DAISIE_MW_loglik_choosepar = function(
  trparsopt,
  trparsfix,
  idparsopt,
  idparsfix,
  pars2,
  datalist,
  methode = 'lsodes',
  CS_version = 1,
  abstolint = 1E-16,
  reltolint = 1E-14,
  distance_type = 'continent',
  distance_dep = 'power',
  parallel = 'local',
  cpus = 3
)
{
  distance_dep_options1 <- c('sigmoidal_col','sigmoidal_ana','sigmoidal_clado','area_additive_clado','area_interactive_clado','area_interactive_clado0','area_interactive_clado1','area_interactive_clado2','area_interactive_clado3')
  trpars1 = rep(0,10 + is.element(distance_dep,distance_dep_options1) + 2 * (distance_dep == 'sigmoidal_col_ana'))
  trpars1[idparsopt] = trparsopt
  if(length(idparsfix) != 0)
  {
    trpars1[idparsfix] = trparsfix
  }
  if(max(trpars1[-5]) >= 1 || trpars1[5] > 1 || min(trpars1) < 0)
  {
    loglik = -Inf
  } else
  {
    pars1 = trpars1/(1 - trpars1)
    if(min(pars1) < 0)
    {
      loglik = -Inf
    } else {
      pars1[c(4,8)] <- -pars1[c(4,8)]
      if(distance_dep == 'sigmoidal_col')
      {
        pars1[8] <- -pars1[8]
      }
      for(i in 1:length(datalist))
      {
        area <- datalist[[i]][[1]]$area
        distance <- (distance_type == 'continent') * datalist[[i]][[1]]$distance_continent +
          (distance_type == 'nearest_big') * datalist[[i]][[1]]$distance_nearest_big +
          (distance_type == 'biologically_realistic') * datalist[[i]][[1]]$distance_biologically_realistic
        if(distance == 0)
        {
          stop('Distance to the mainland is 0 in the data.')
        }
        if(distance_dep == 'exp')
        {
          distance <- exp(distance)
        }
        pars1new = pars1[c(1,3,5,7,9)] * c(rep(area,3),rep(distance,2))^pars1[c(2,4,6,8,10)]
        if(distance_dep == 'sigmoidal_col')
        {
          pars1new[4] <- f_sigmoidal_neg(d = distance,k = pars1[7],x = pars1[8],d0 = pars1[11])
        } else
        if(distance_dep == 'sigmoidal_ana')
        {
          pars1new[5] <- f_sigmoidal_pos(d = distance,k = pars1[9],x = pars1[10],d0 = pars1[11])
        } else
        if(distance_dep == 'sigmoidal_clado')
        {
          pars1new[1] <- f_sigmoidal_pos(d = distance,k = pars1[1],x = pars1[2],d0 = pars1[11])
        } else
        if(distance_dep == 'sigmoidal_col_ana')
        {
          pars1new[4] <- f_sigmoidal_neg(d = distance,k = pars1[7],x = pars1[8],d0 = pars1[11])
          pars1new[5] <- f_sigmoidal_pos(d = distance,k = pars1[9],x = pars1[10],d0 = pars1[12])
        } else
        if(distance_dep == 'area_additive_clado')
        {
          pars1new[1] <- pars1[1] * area^pars1[2] * distance^pars1[11]
        } else
        if(distance_dep == 'area_interactive_clado' || distance_dep == 'area_interactive_clado0')
        {
          pars1new[1] <- pars1[1] * area^(pars1[2] + pars1[11] * log(distance))
        } else
        if(distance_dep == 'area_interactive_clado1')
        {
          pars1new[1] <- pars1[1] * area^(pars1[2] + distance/pars1[11])
        } else
        if(distance_dep == 'area_interactive_clado2')
        {
          pars1new[1] <- pars1[1] * area^(pars1[2] + 1/(1 + pars1[11]/distance))
        } else
        if(distance_dep == 'area_interactive_clado3')
        {
          pars1new[1] <- pars1[1] * (area + distance/pars1[11])^pars1[2]
        } else
        if(distance_dep != 'power' && distance_dep != 'exp')
        {
          stop('This type of distance_dep is not supported. Please check spelling.')
        }
        M = datalist[[i]][[1]]$not_present + length(datalist[[i]]) - 1
        pars1new[4] = pars1new[4] / M
        datalist[[i]][[1]]$pars1new = pars1new
      }
      if(max(pars1new[c(1,2,4,5)]) > 1E+10)
      {
        loglik = -Inf
      } else
      {
        if(parallel != 'no')
        {
          if(parallel == 'local')
          {
            cpus = min(cpus,parallel::detectCores())
            cl = parallel::makeCluster(cpus - 1)
            doParallel::registerDoParallel(cl)
            on.exit(parallel::stopCluster(cl))
          } else
            if(parallel == 'cluster')
            {
              if(.Platform$OS.type != "unix")
              {
                cat('cluster does not work on a non-unix environment, choose local instead.\n')
                return(-Inf)
              }
              doMC::registerDoMC(cpus - 1)
            }
          X = NULL; rm(X)
          loglik = foreach::foreach(X = datalist,.combine = sum,.export = c("pars2"),.packages = c('DAISIE','foreach','desolve','doParallel'))  %dopar%  DAISIE_loglik_all(X[[1]]$pars1new,pars2,X)
        } else {
          loglik = 0
          if(pars2[4] == 0.5) pb <- utils::txtProgressBar(min = 0, max = length(datalist), style = 3)
          for(i in 1:length(datalist))
          {
            loglik = loglik + DAISIE_loglik_all(pars1 = datalist[[i]][[1]]$pars1new,pars2 = pars2,datalist = datalist[[i]],methode = methode,CS_version,abstolint = abstolint,reltolint = reltolint)
            if(pars2[4] == 0.5) utils::setTxtProgressBar(pb, i)
          }
          if(pars2[4] == 0.5) close(pb)
        }
      }
      if(is.nan(loglik) || is.na(loglik))
      {
        cat("There are parameter values used which cause numerical problems.\n")
        loglik = -Inf
      }
    }
  }
  return(loglik)
}


#' @name DAISIE_MW_ML
#' @title Maximization of the loglikelihood under the DAISIE model with clade-specific
#' diversity-dependence and explicit dependencies on island area and distance
#' from the mainland or nearest landmass as hypothesized by MacArthur & Wilson
#' @description This function computes the maximum likelihood estimates of the parameters of
#' the relationships between parameters of the DAISIE model with clade-specific
#' diversity-dependence and island area and distance of the island to the
#' mainlandor nearest landmass, for data from lineages colonizing several
#' islands/archipelagos. It also outputs the corresponding loglikelihood that
#' can be used in model comparisons.
#' 
#' A note on the sigmoidal functions used in distance_dep: For anagenesis and
#' cladogenesis, the functional relationship is k * (d/d0)^x/(1 + (d/d0)^x);
#' for colonization the relationship is: k - k * (d/d0)^x/(1 + (d/d0)^x). The
#' d0 parameter is the 11th parameter entered. In the of 'sigmoidal_col_ana',
#' the 11th parameter is the d0 for colonization and the 12th is the d0 for
#' anagenesis.
#' 
#' @param datalist Data object containing information on colonisation and
#' branching times. This object can be generated using the DAISIE_dataprep
#' function, which converts a user-specified data table into a data object, but
#' the object can of course also be entered directly. It is an R list object
#' with the following elements.\cr The first element of the list has two three
#' components: \cr \cr \code{$island_age} - the island age \cr Then, depending
#' on whether a distinction between types is made, we have:\cr
#' \code{$not_present} - the number of mainland lineages that are not present
#' on the island \cr\cr The remaining elements of the list each contains
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
#' applicable for endemic clades)
#' @param initparsopt The initial values of the parameters that must be
#' optimized; they are all positive
#' @param idparsopt The ids of the parameters that must be optimized. The ids
#' are defined as follows: \cr \cr id = 1 corresponds to lambda^c0
#' (cladogenesis rate for unit area) \cr id = 2 corresponds to y (exponent of
#' area for cladogenesis rate) \cr id = 3 corresponds to mu0 (extinction rate
#' for unit area) \cr id = 4 corresponds to x (exponent of 1/area for
#' extinction rate) \cr id = 5 corresponds to K0 (clade-level carrying capacity
#' for unit area) \cr id = 6 corresponds to z (exponent of area for clade-level
#' carrying capacity) \cr id = 7 corresponds to gamma0 (immigration rate for
#' unit distance) \cr id = 8 corresponds to alpha (exponent of 1/distance for
#' immigration rate) \cr id = 9 corresponds to lambda^a0 (anagenesis rate for
#' unit distance) \cr id = 10 corresponds to beta (exponent of 1/distance for
#' anagenesis rate) \cr
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
#' @param distance_type Use 'continent' if the distance to the continent should
#' be used, use 'nearest_big' if the distance to the nearest big landmass
#' should be used, and use 'biologically_realistic' if the distance should take
#' into account some biologically realism, e.g. an average of the previous two
#' if both are thought to contribute.
#' @param distance_dep Sets what type of distance dependence should be used.
#' Default is a power law, denoted as 'power'. Alternatives are an exponantial
#' relationship denoted by 'exp' or sigmoids, either 'sigmoidal_col' for a
#' sigmoid in the colonization, 'sigmoidal_ana' for sigmoidal anagenesis,
#' 'sigmoidal_clado' for sigmoidal cladogenesis, and 'sigmoidal_col_ana' for
#' signoids in both colonization and anagenesis.
#' @param parallel Sets whether parallel computation should be used. Use 'no'
#' if no parallel computing should be used, 'cluster' for parallel computing on
#' a unix/linux cluster, and 'local' for parallel computation on a local
#' machine.
#' @param cpus Number of cpus used in parallel computing. Default is 3. Will
#' not have an effect if parallel = 'no'.
#' @return The output is a dataframe containing estimated parameters and
#' maximum loglikelihood.  \item{lambda_c0}{ gives the maximum likelihood
#' estimate of lambda^c, the rate of cladogenesis for unit area} \item{y}{
#' gives the maximum likelihood estimate of y, the exponent of area for the
#' rate of cladogenesis} \item{mu0}{ gives the maximum likelihood estimate of
#' mu0, the extinction rate} \item{x}{ gives the maximum likelihood estimate of
#' x, the exponent of 1/area for the extinction rate} \item{K0}{ gives the
#' maximum likelihood estimate of K0, the carrying-capacity for unit area}
#' \item{z}{ gives the maximum likelihood estimate of z, the exponent of area
#' for the carrying capacity} \item{gamma0}{ gives the maximum likelihood
#' estimate of gamma0, the immigration rate for unit distance} \item{y}{ gives
#' the maximum likelihood estimate of alpha, the exponent of 1/distance for the
#' rate of colonization} \item{lambda_a0}{ gives the maximum likelihood
#' estimate of lambda^a0, the rate of anagenesis for unit distance}
#' \item{beta}{ gives the maximum likelihood estimate of beta, the exponent of
#' 1/distance for the rate of anagenesis} \item{loglik}{ gives the maximum
#' loglikelihood} \item{df}{ gives the number of estimated parameters, i.e.
#' degrees of feedom} \item{conv}{ gives a message on convergence of
#' optimization; conv = 0 means convergence}
#' @author Rampal S. Etienne
#' @seealso \code{\link{DAISIE_ML_CS}},
#' @references Valente, L.M., A.B. Phillimore and R.S. Etienne (2015).
#' Equilibrium and non-equilibrium dynamics simultaneously operate in the
#' Galapagos islands. Ecology Letters 18: 844-852. <DOI:10.1111/ele.12461>.
#' @keywords models
#' @examples
#' 
#' cat('No examples')
#' 
#' @export DAISIE_MW_ML
DAISIE_MW_ML = function(
  datalist,
  initparsopt,
  idparsopt,
  parsfix,
  idparsfix,
  res = 100,
  ddmodel = 11,
  cond = 0,
  island_ontogeny = NA,
  tol = c(1E-4, 1E-5, 1E-7),
  maxiter = 1000 * round((1.25)^length(idparsopt)),
  methode = "lsodes",
  optimmethod = 'subplex',
  CS_version = 1,
  verbose = 0,
  tolint = c(1E-16,1E-10),
  distance_type = 'continent',
  distance_dep = 'power',
  parallel = 'local',
  cpus = 3
)
{
  options(warn=-1)
  distance_dep_options1 <- c('sigmoidal_col','sigmoidal_ana','sigmoidal_clado','area_additive_clado','area_interactive_clado','area_interactive_clado0','area_interactive_clado1','area_interactive_clado2','area_interactive_clado3')
  numpars <- 10 + is.element(distance_dep,distance_dep_options1) + 2 * (distance_dep == 'sigmoidal_col_ana')
  if(numpars == 11)
  {
    out2err = data.frame(lambda_c0 = NA, y = NA, mu_0 = NA, x = NA, K_0 = NA, z = NA, gamma_0 = NA, alpha = NA, lambda_a0 = NA, beta = NA, d0 = NA, loglik = NA, df = NA, conv = NA)
  } else
  if(numpars == 12)
  {
    out2err = data.frame(lambda_c0 = NA, y = NA, mu_0 = NA, x = NA, K_0 = NA, z = NA, gamma_0 = NA, alpha = NA, lambda_a0 = NA, beta = NA, d0_col = NA, d0_ana = NA, loglik = NA, df = NA, conv = NA)
  } else
  {
    out2err = data.frame(lambda_c0 = NA, y = NA, mu_0 = NA, x = NA, K_0 = NA, z = NA, gamma_0 = NA, alpha = NA, lambda_a0 = NA, beta = NA, loglik = NA, df = NA, conv = NA)
  }
  out2err = invisible(out2err)
  idpars = sort(c(idparsopt,idparsfix))
  if((prod(idpars == (1:numpars)) != 1) || (length(initparsopt) != length(idparsopt)) || (length(parsfix) != length(idparsfix)))
  {
    cat("The parameters to be optimized and/or fixed are incoherent.\n")
    return(out2err)
  }
  if(length(idparsopt) > numpars)
  {
    cat("The number of parameters to be optimized is too high.\n")
    return(out2err)
  }
  namepars = c("lambda_c0","y","mu_0","x","K_0","z","gamma_0","alpha","lambda_a0","beta")
  if(is.element(distance_dep,distance_dep_options1))
  {
    namepars = c(namepars,"d0")   
  } else if(distance_dep == 'sigmoidal_col_ana')
  {
    namepars = c(namepars,"d0_col","d0_ana")   
  }
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
  pars2 = c(res,ddmodel,cond,verbose,island_ontogeny)
  optimpars = c(tol,maxiter)
  initloglik = DAISIE_MW_loglik_choosepar(trparsopt = trparsopt,trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,pars2 = pars2,datalist = datalist,methode = methode,CS_version = CS_version,abstolint = tolint[1],reltolint = tolint[2],distance_type = distance_type,parallel = parallel,cpus = cpus,distance_dep = distance_dep)
  cat("The loglikelihood for the initial parameter values is",initloglik,"\n")
  if(initloglik == -Inf)
  {
    cat("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n")
    return(out2err)
  }
  cat("Optimizing the likelihood - this may take a while.","\n")
  utils::flush.console()
  out = DDD::optimizer(optimmethod = optimmethod,optimpars = optimpars,fun = DAISIE_MW_loglik_choosepar,trparsopt = trparsopt,idparsopt = idparsopt,trparsfix = trparsfix,idparsfix = idparsfix,pars2 = pars2,datalist = datalist,methode = methode,CS_version = CS_version,abstolint = tolint[1],reltolint = tolint[2],distance_type = distance_type,parallel = parallel,cpus = cpus,distance_dep = distance_dep)
  if(out$conv != 0)
  {
    cat("Optimization has not converged. Try again with different initial values.\n")
    out2 = out2err
    out2$conv = out$conv
    return(out2)
  }
  MLtrpars = as.numeric(unlist(out$par))
  MLpars = MLtrpars/(1 - MLtrpars)
  ML = as.numeric(unlist(out$fvalues))
  MLpars1 = rep(0,numpars)
  MLpars1[idparsopt] = MLpars
  if(length(idparsfix) != 0) { MLpars1[idparsfix] = parsfix }
  if(MLpars1[5] > 10^7){ MLpars1[5] = Inf }
  s1output <- function(MLpars1,distance_dep)
  {
     s1 <- switch(distance_dep,
     power = sprintf('Maximum likelihood parameter estimates:\n
               lambda_c = %f * A^ %f\n
               mu = %f * A^ -%f\n
               K = %f * A^ %f\n
               M * gamma = %f * d^ -%f\n
               lambda_a = %f * d^ %f\n',MLpars1[1],MLpars1[2],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6],MLpars1[7],MLpars1[8],MLpars1[9],MLpars1[10]),
     signoidal_col = sprintf('Maximum likelihood parameter estimates:\n
               lambda_c = %f * A^ %f\n
               mu = %f * A^ -%f\n
               K = %f * A^ %f\n
               M * gamma = %f * (1 - (d/%f)^%f / (1 + (d/%f)^%f )\n
               lambda_a = %f * d^ %f\n',MLpars1[1],MLpars1[2],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6],MLpars1[7],MLpars1[11],MLpars1[8],MLpars1[11],MLpars1[8],MLpars1[9],MLpars1[10]),
     sigmoidal_ana = sprintf('Maximum likelihood parameter estimates:\n
               lambda_c = %f * A^ %f\n
               mu = %f * A^ -%f\n
               K = %f * A^ %f\n
               M * gamma = %f * d^ -%f\n
               lambda_a = %f * (d/%f)^%f / (1 + (d/%f)^%f )\n',MLpars1[1],MLpars1[2],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6],MLpars1[7],MLpars1[8],MLpars1[9],MLpars1[11],MLpars1[10],MLpars1[11],MLpars1[10]),
     sigmoidal_clado = sprintf('Maximum likelihood parameter estimates:\n
               lambda_c = %f * (d/%f)^%f / (1 + (d/%f)^%f )\n
               mu = %f * A^ -%f\n
               K = %f * A^ %f\n
               M * gamma = %f * d^ -%f\n
               lambda_a = %f * d^ %f\n',MLpars1[1],MLpars1[11],MLpars1[2],MLpars1[11],MLpars1[2],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6],MLpars1[7],MLpars1[8],MLpars1[9],MLpars1[10]),
     sigmoidal_col_ana = sprintf('Maximum likelihood parameter estimates:\n
               lambda_c = %f * A^ %f\n
               mu = %f * A^ -%f\n
               K = %f * A^ %f\n
               M * gamma = %f * (1 - (d/%f)^%f / (1 + (d/%f)^%f )\n
               lambda_a = %f * (d/%f)^%f / (1 + (d/%f)^%f )\n',MLpars1[1],MLpars1[2],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6],MLpars1[7],MLpars1[11],MLpars1[8],MLpars1[11],MLpars1[8],MLpars1[9],MLpars1[12],MLpars1[10],MLpars1[12],MLpars1[10]),
     area_additive_clado = sprintf('Maximum likelihood parameter estimates:\n
               lambda_c = %f * A^ %f * d^ %f \n
               mu = %f * A^ -%f\n
               K = %f * A^ %f\n
               M * gamma = %f * d^ -%f\n
               lambda_a = %f * d^ %f\n',MLpars1[1],MLpars1[2],MLpars1[11],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6],MLpars1[7],MLpars1[8],MLpars1[9],MLpars1[10]),
     area_interactive_clado = sprintf('Maximum likelihood parameter estimates:\n
               lambda_c = %f * A^ (%f + %f * log(d)) \n
               mu = %f * A^ -%f\n
               K = %f * A^ %f\n
               M * gamma = %f * d^ -%f\n
               lambda_a = %f * d^ %f\n',MLpars1[1],MLpars1[2],MLpars1[11],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6],MLpars1[7],MLpars1[8],MLpars1[9],MLpars1[10]),
     area_interactive_clado0 = sprintf('Maximum likelihood parameter estimates:\n
               lambda_c = %f * A^ (%f + %f * log(d)) \n
               mu = %f * A^ -%f\n
               K = %f * A^ %f\n
               M * gamma = %f * d^ -%f\n
               lambda_a = %f * d^ %f\n',MLpars1[1],MLpars1[2],MLpars1[11],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6],MLpars1[7],MLpars1[8],MLpars1[9],MLpars1[10]),    
     area_interactive_clado1 = sprintf('Maximum likelihood parameter estimates:\n
               lambda_c = %f * A^ (%f + d/%f)) \n
               mu = %f * A^ -%f\n
               K = %f * A^ %f\n
               M * gamma = %f * d^ -%f\n
               lambda_a = %f * d^ %f\n',MLpars1[1],MLpars1[2],MLpars1[11],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6],MLpars1[7],MLpars1[8],MLpars1[9],MLpars1[10]),
     area_interactive_clado2 = sprintf('Maximum likelihood parameter estimates:\n
               lambda_c = %f * A^ (%f + d/(d + %f)) \n
               mu = %f * A^ -%f\n
               K = %f * A^ %f\n
               M * gamma = %f * d^ -%f\n
               lambda_a = %f * d^ %f\n',MLpars1[1],MLpars1[2],MLpars1[11],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6],MLpars1[7],MLpars1[8],MLpars1[9],MLpars1[10]),
     area_interactive_clado3 = sprintf('Maximum likelihood parameter estimates:\n
               lambda_c = %f * (A + d/%f)^ %f\n \n
               mu = %f * A^ -%f\n
               K = %f * A^ %f\n
               M * gamma = %f * d^ -%f\n
               lambda_a = %f * d^ %f\n',MLpars1[1],MLpars1[11],MLpars1[2],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6],MLpars1[7],MLpars1[8],MLpars1[9],MLpars1[10])
     )
     return(s1)
  }
  s2 = sprintf('Maximum loglikelihood: %f',ML)
  if(is.element(distance_dep,distance_dep_options1))
  {
    out2 = data.frame(lambda_c0 = MLpars1[1], y = MLpars1[2], mu_0 = MLpars1[3], x = MLpars1[4], K_0 = MLpars1[5], z = MLpars1[6], gamma_0 = MLpars1[7], alpha = MLpars1[8], lambda_a0 = MLpars1[9], beta = MLpars1[10], d_0 = MLpars1[11], loglik = ML, df = length(initparsopt), conv = unlist(out$conv))
  } else
  if(distance_dep == 'sigmoidal_col_ana')
  {
    out2 = data.frame(lambda_c0 = MLpars1[1], y = MLpars1[2], mu_0 = MLpars1[3], x = MLpars1[4], K_0 = MLpars1[5], z = MLpars1[6], gamma_0 = MLpars1[7], alpha = MLpars1[8], lambda_a0 = MLpars1[9], beta = MLpars1[10], d0_col = MLpars1[11], d0_ana = MLpars1[12], loglik = ML, df = length(initparsopt), conv = unlist(out$conv))
  } else
  {
    out2 = data.frame(lambda_c0 = MLpars1[1], y = MLpars1[2], mu_0 = MLpars1[3], x = MLpars1[4], K_0 = MLpars1[5], z = MLpars1[6], gamma_0 = MLpars1[7], alpha = MLpars1[8], lambda_a0 = MLpars1[9], beta = MLpars1[10], loglik = ML, df = length(initparsopt), conv = unlist(out$conv))
  }
  cat("\n",s1output(MLpars1,distance_dep),"\n",s2,"\n")
  return(invisible(out2))
}
