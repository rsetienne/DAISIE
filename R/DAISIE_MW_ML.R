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

convert_parameters_MW <- function(pars1,area,distance,M,distance_dep) {
  pars1new = pars1[c(1,3,5,7,9)] * c(rep(area,3),rep(distance,2))^pars1[c(2,4,6,8,10)]
  if(distance_dep == 'sigmoidal_col')
  {
    pars1new[4] <- f_sigmoidal_neg(d = distance,k = pars1[7],x = pars1[8],d0 = pars1[11])
  } else
    if(distance_dep == 'sigmoidal_ana')
    {
      pars1new[5] <- f_sigmoidal_pos(d = distance,k = pars1[9],x = pars1[10],d0 = pars1[11])
    } else if(distance_dep == 'sigmoidal_clado')
    {
      pars1new[1] <- f_sigmoidal_pos(d = distance,k = pars1[1],x = pars1[2],d0 = pars1[11])
    } else if(distance_dep == 'sigmoidal_col_ana')
    {
      pars1new[4] <- f_sigmoidal_neg(d = distance,k = pars1[7],x = pars1[8],d0 = pars1[11])
      pars1new[5] <- f_sigmoidal_pos(d = distance,k = pars1[9],x = pars1[10],d0 = pars1[12])
    } else if(distance_dep == 'area_additive_clado')
    {
      pars1new[1] <- pars1[1] * area^pars1[2] * distance^pars1[11]
    } else if(distance_dep == 'area_interactive_clado' || distance_dep == 'area_interactive_clado0')
    {
      pars1new[1] <- pars1[1] * area^(pars1[2] + pars1[11] * log(distance))
    } else if(distance_dep == 'area_interactive_clado1')
    {
      pars1new[1] <- pars1[1] * area^(pars1[2] + distance/pars1[11])
    } else if(distance_dep == 'area_interactive_clado2')
    {
      pars1new[1] <- pars1[1] * area^(pars1[2] + 1/(1 + pars1[11]/distance))
    } else if(distance_dep == 'area_interactive_clado3')
    {
      pars1new[1] <- pars1[1] * (area + distance/pars1[11])^pars1[2]
    } else if(distance_dep == 'area_interactive_clado4')
    {
      pars1new[1] <- pars1[1] * area^(pars1[2]/(1 + pars1[11]/distance))
    } else if(distance_dep != 'power' && distance_dep != 'exp')
    {
      stop('This type of distance_dep is not supported. Please check spelling.')
    }
  pars1new[4] = pars1new[4] / M
  return(pars1new)
}

distance_dep_options1_fun <- function()
{
  return(c('sigmoidal_col',
           'sigmoidal_ana',
           'sigmoidal_clado',
           'area_additive_clado',
           'area_interactive_clado',
           'area_interactive_clado0',
           'area_interactive_clado1',
           'area_interactive_clado2',
           'area_interactive_clado3',
           'area_interactive_clado4')
  )
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
  methode = 'odeint::runge_kutta_cash_karp54',
  CS_version = list(model = 1, function_to_optimize = 'DAISIE'),
  abstolint = 1E-16,
  reltolint = 1E-14,
  distance_type = 'continent',
  distance_dep = 'power',
  parallel = 'local',
  cpus = 3
)
{
  distance_dep_options1 <- distance_dep_options1_fun()
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
        if(distance_type == 'simulation')
        {
          distance <- datalist[[i]][[1]]$distance
        } else {
          distance <- (distance_type == 'continent') * datalist[[i]][[1]]$distance_continent +
            (distance_type == 'nearest_big') * datalist[[i]][[1]]$distance_nearest_big +
            (distance_type == 'biologically_realistic') * datalist[[i]][[1]]$distance_biologically_realistic
        }
        if(distance == 0)
        {
          stop('Distance to the mainland is 0 in the data.')
        }
        if(distance_dep == 'exp')
        {
          distance <- exp(distance)
        }
        M = datalist[[i]][[1]]$not_present + length(datalist[[i]]) - 1
        pars1new <- convert_parameters_MW(pars1 = pars1,area = area,distance = distance,M = M,distance_dep = distance_dep)
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
                warning('cluster does not work on a non-unix environment, choose local instead.')
                return(-Inf)
              }
              doMC::registerDoMC(cpus - 1)
            }
          X = NULL; rm(X)

          suppressWarnings({
            loglik = foreach::foreach(
              X = datalist,
              .combine = sum,
              .export = c("pars2"),
              .packages = c('DAISIE','foreach','deSolve','doParallel')) %dopar%
              DAISIE_loglik_all(pars1 = X[[1]]$pars1new,
                                pars2 = pars2,
                                datalist = X,
                                methode = methode,
                                CS_version = CS_version,
                                abstolint = abstolint,
                                reltolint = reltolint)
          })
        } else {
          loglik = 0
          if(pars2[4] == 0.5) pb <- utils::txtProgressBar(min = 0, max = length(datalist), style = 3)
          for(i in 1:length(datalist))
          {
            loglik = loglik + DAISIE_loglik_all(pars1 = datalist[[i]][[1]]$pars1new,
                                                pars2 = pars2,datalist = datalist[[i]],
                                                methode = methode,
                                                CS_version = CS_version,
                                                abstolint = abstolint,
                                                reltolint = reltolint)
            if(pars2[4] == 0.5) utils::setTxtProgressBar(pb, i)
          }
          if(pars2[4] == 0.5) close(pb)
        }
      }
      if(is.nan(loglik) || is.na(loglik))
      {
        warning("There are parameter values used which cause numerical problems.")
        loglik = -Inf
      }
    }
  }
  return(loglik)
}


#' @name DAISIE_MW_ML
#' @title Maximization of the loglikelihood under the DAISIE model with clade-specific
#' diversity-dependence and explicit dependencies on island area and isolation
#' as hypothesized by MacArthur & Wilson
#' @description This function computes the maximum likelihood estimates of the parameters of
#' the relationships between parameters of the DAISIE model (with clade-specific
#' diversity-dependence) and island area and distance of the island to the
#' mainland for data from lineages colonizing several
#' islands/archipelagos. It also outputs the corresponding loglikelihood that
#' can be used in model comparisons.
#'
#' A note on the sigmoidal functions used in distance_dep: For anagenesis and
#' cladogenesis, the functional relationship is k * (d/d0)^x/(1 + (d/d0)^x);
#' for colonization the relationship is: k - k * (d/d0)^x/(1 + (d/d0)^x). The
#' d0 parameter is the 11th parameter entered. In 'sigmoidal_col_ana',
#' the 11th parameter is the d0 for colonization and the 12th is the d0 for
#' anagenesis.
#'
#' @inheritParams default_params_doc
#' @param datalist Data object containing information on colonisation and
#' branching times of species for several islands or archipelagos, as well as the area,
#' isolation and age of each of the islands/archipelagos. See data(archipelagos41) for
#' an example.
#' @param initparsopt The initial values of the parameters that must be
#' optimized; they are all positive
#' @param idparsopt The ids of the parameters that must be optimized. The ids
#' are defined as follows (see Valente et al 2020 Supplementary Tables 1 and 2
#' a better explanation of the models and parameters): \cr \cr
#' id = 1 corresponds to lambda^c0 (cladogenesis rate for unit area) \cr
#' id = 2 corresponds to y (exponent of area for cladogenesis rate) \cr
#' id = 3 corresponds to mu0 (extinction rate for unit area) \cr
#' id = 4 corresponds to x (exponent of 1/area for extinction rate) \cr
#' id = 5 corresponds to K0 (clade-level carrying capacity for unit area) \cr
#' id = 6 corresponds to z (exponent of area for clade-level carrying capacity) \cr
#' id = 7 corresponds to gamma0 (immigration rate for unit distance) \cr
#' id = 8 corresponds to alpha (exponent of 1/distance for immigration rate) \cr id = 9 corresponds to lambda^a0 (anagenesis rate for
#' unit distance) \cr
#' id = 10 corresponds to beta (exponent of 1/distance for anagenesis rate) \cr
#' id = 11 corresponds to d0 in models M15 to M19, and models with
#' distance_dep = 'sigmoidal_col', 'sigmoidal_ana' or 'sigmoidal_clado';
#' or d0 for colonisation (when specifying distance_dep = 'sigmoidal_col_ana'\cr
#' id = 12 corresponds to d0 for anagenesis when specifying
#' distance_dep = 'sigmoidal_col_ana' \cr
#' @param idparsfix The ids of the parameters that should not be optimized,
#' e.g. c(1,3) if lambda^c and K should not be optimized.
#' @param parsfix The values of the parameters that should not be optimized
#' @param res Sets the maximum number of species for which a probability must
#' be computed, must be larger than the size of the largest clade
#' @param ddmodel Sets the model of diversity-dependence: \cr \cr ddmodel = 0 :
#' no diversity dependence \cr
#' ddmodel = 1 : linear dependence in speciation rate \cr
#' ddmodel = 11: linear dependence in speciation rate and in immigration rate \cr
#' ddmodel = 2 : exponential dependence in speciation rate\cr
#' ddmodel = 21: exponential dependence in speciation rate and in immigration rate\cr
#' @param cond cond = 0 : conditioning on island age \cr
#' cond = 1 : conditioning on island age and non-extinction of the island
#' biota \cr
#' @param island_ontogeny type of island ontonogeny. If NA, then constant ontogeny is assumed
#' @param tol Sets the tolerances in the optimization. Consists of: \cr
#' reltolx = relative tolerance of parameter values in optimization \cr
#' reltolf = relative tolerance of function value in optimization \cr
#' abstolx = absolute tolerance of parameter values in optimization
#' @param maxiter Sets the maximum number of iterations in the optimization
#' @param methode Method of the ODE-solver. See package deSolve for details.
#' Default is "lsodes"
#' @param optimmethod Method used in likelihood optimization. Default is
#' "subplex" (see subplex package). Alternative is 'simplex' which was the
#' method in previous versions.
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
#' Default is a power law, denoted as 'power' (models M1-14 in Valente et al 2020).
#' Alternatives are additive or interactive contributions of distance and area
#' to the rate of cladogenesis ("area_additive_clado"; "area_interactive_clado",
#' "area_interactive_clado1" and "area_interactive_clado2"). Other alternatives are
#'  exponential relationship denoted by 'exp'; or sigmoids, either
#' 'sigmoidal_col' for a sigmoid in the colonization, 'sigmoidal_ana' for sigmoidal anagenesis,
#' 'sigmoidal_clado' for sigmoidal cladogenesis, and 'sigmoidal_col_ana' for
#' sigmoids in both colonization and anagenesis. \cr
#' A key for the different options of distance_dep that should be specified to run the
#' models from Valente et al 2020 (Supplementary Data Table 1 and 2) is given below: \cr
#' * M1 to M14 - 'power' \cr
#' * M15 -'area_additive_clado' \cr
#' * M16 and M19 -'area_interactive_clado' \cr
#' * M17 -'area_interactive_clado1' \cr
#' * M18 - 'area_interactive_clado2' \cr
#' * M20 and M24 - sigmoidal_col' \cr
#' * M21, M25 and M28 - sigmoidal_ana' \cr
#' * M22 and M26 - 'sigmoidal_clado'\cr
#' * M23 and M27 - 'sigmoidal_col_ana' \cr
#' @param parallel Sets whether parallel computation should be used. Use 'no'
#' if no parallel computing should be used, 'cluster' for parallel computing on
#' a unix/linux cluster, and 'local' for parallel computation on a local
#' machine.
#' @param cpus Number of cpus used in parallel computing. Default is 3. Will
#' not have an effect if parallel = 'no'.
#' @param num_cycles The number of cycles the optimizer will go through.
#'   Default is 1.
#' @return The output is a dataframe containing estimated parameters and
#' maximum loglikelihood.
#' \item{lambda_c0}{ gives the maximum likelihood estimate of lambda^c,
#' the rate of cladogenesis for unit area}
#' \item{y}{ gives the maximum likelihood estimate of y,
#' the exponent of area for the rate of cladogenesis}
#' \item{mu0}{ gives the maximum likelihood estimate of mu0,
#' the extinction rate}
#' \item{x}{ gives the maximum likelihood estimate of
#' x, the exponent of 1/area for the extinction rate}
#' \item{K0}{ gives the maximum likelihood estimate of K0,
#' the carrying-capacity for unit area}
#' \item{z}{ gives the maximum likelihood estimate of z, the exponent of area
#' for the carrying capacity}
#' \item{gamma0}{ gives the maximum likelihood
#' estimate of gamma0, the immigration rate for unit distance} \item{y}{ gives
#' the maximum likelihood estimate of alpha, the exponent of 1/distance for the
#' rate of colonization}
#' \item{lambda_a0}{ gives the maximum likelihood
#' estimate of lambda^a0, the rate of anagenesis for unit distance}
#' \item{beta}{ gives the maximum likelihood estimate of beta, the exponent of
#' 1/distance for the rate of anagenesis}
#' \item{loglik}{ gives the maximum loglikelihood}
#' \item{df}{ gives the number of estimated parameters, i.e. degrees of feedom}
#' \item{conv}{ gives a message on convergence of optimization;
#' conv = 0 means convergence}
#' @author Rampal S. Etienne & Luis Valente
#' @seealso \code{\link{DAISIE_ML_CS}},
#' @references Valente L, Phillimore AB, Melo M, Warren BH, Clegg SM, Havenstein K,
#'  Tiedemann R, Illera JC, Thébaud C, Aschenbach T, Etienne RS. A simple dynamic model
#'  explains island bird diversity worldwide (2020) Nature, 579, 92-96
#' @keywords models
#' @examples
#'
#' cat("
#' ### Fit the M19 model as in Valente et al 2020, using the ML
#' parameters as starting values (see Supplementary Tables 1 and 2).
#'
#' utils::data(archipelagos41)
#'
#' DAISIE_MW_ML(
#' datalist= archipelagos41,
#' initparsopt =
#' c(0.040073803,	1.945656546,	0.150429656,
#' 67.25643672,	0.293635061,	0.059096872,	0.382688527,
#' 0.026510781),
#' idparsopt = c(1,3,4,7,8,9,10,11),
#' parsfix = c(0,Inf,0) ,
#' idparsfix = c(2,5,6),
#' res = 100,
#' ddmodel = 0,
#' methode = 'lsodes',
#' cpus = 4,
#' parallel = 'local',
#' optimmethod = 'subplex',
#' tol = c(1E-4, 1E-5, 1E-7),
#' distance_type = 'continent',
#' distance_dep = 'area_interactive_clado'
#' )
#'")
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
  methode = "odeint::runge_kutta_cash_karp54",
  optimmethod = 'simplex',
  CS_version = list(model = 1, function_to_optimize = 'DAISIE'),
  verbose = 0,
  tolint = c(1E-16, 1E-10),
  distance_type = 'continent',
  distance_dep = 'power',
  parallel = 'local',
  cpus = 3,
  num_cycles = 1
)
{
  distance_dep_options1 <- distance_dep_options1_fun()
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
    warning("The parameters to be optimized and/or fixed are incoherent.\n")
    return(out2err)
  }
  if(length(idparsopt) > numpars)
  {
    warning("The number of parameters to be optimized is too high.\n")
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

  print_ml_par_settings(
    namepars = namepars,
    idparsopt = idparsopt,
    idparsfix = idparsfix,
    idparsnoshift = NA,
    all_no_shift = NA,
    verbose = verbose
  )

  trparsopt = initparsopt/(1 + initparsopt)
  trparsopt[which(initparsopt == Inf)] = 1
  trparsfix = parsfix/(1 + parsfix)
  trparsfix[which(parsfix == Inf)] = 1
  pars2 = c(res,ddmodel,cond,verbose,island_ontogeny)
  optimpars = c(tol,maxiter)
  initloglik = DAISIE_MW_loglik_choosepar(trparsopt = trparsopt,trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,pars2 = pars2,datalist = datalist,methode = methode,CS_version = CS_version,abstolint = tolint[1],reltolint = tolint[2],distance_type = distance_type,parallel = parallel,cpus = cpus,distance_dep = distance_dep)
  message("The loglikelihood for the initial parameter values is ", initloglik)
  if(initloglik == -Inf)
  {
    warning("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n")
    return(out2err)
  }
  message("Optimizing the likelihood - this may take a while.")
  out = DDD::optimizer(optimmethod = optimmethod,
                       optimpars = optimpars,
                       fun = DAISIE_MW_loglik_choosepar,
                       trparsopt = trparsopt,
                       idparsopt = idparsopt,
                       trparsfix = trparsfix,
                       idparsfix = idparsfix,
                       pars2 = pars2,
                       datalist = datalist,
                       methode = methode,
                       CS_version = CS_version,
                       abstolint = tolint[1],
                       reltolint = tolint[2],
                       distance_type = distance_type,
                       parallel = parallel,
                       cpus = cpus,
                       distance_dep = distance_dep,
                       num_cycles = num_cycles)
  if(out$conv != 0)
  {
    warning("Optimization has not converged. Try again with different initial values.\n")
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
  if(is.element(distance_dep,distance_dep_options1))
  {
    out2 = data.frame(lambda_c0 = MLpars1[1], y = MLpars1[2], mu_0 = MLpars1[3], x = MLpars1[4], K_0 = MLpars1[5], z = MLpars1[6], gamma_0 = MLpars1[7], alpha = MLpars1[8], lambda_a0 = MLpars1[9], beta = MLpars1[10], d_0 = MLpars1[11], loglik = ML, df = length(initparsopt), conv = unlist(out$conv))
  } else {
    if(distance_dep == 'sigmoidal_col_ana')
    {
      out2 = data.frame(lambda_c0 = MLpars1[1], y = MLpars1[2], mu_0 = MLpars1[3], x = MLpars1[4], K_0 = MLpars1[5], z = MLpars1[6], gamma_0 = MLpars1[7], alpha = MLpars1[8], lambda_a0 = MLpars1[9], beta = MLpars1[10], d0_col = MLpars1[11], d0_ana = MLpars1[12], loglik = ML, df = length(initparsopt), conv = unlist(out$conv))
    } else
    {
      out2 = data.frame(lambda_c0 = MLpars1[1], y = MLpars1[2], mu_0 = MLpars1[3], x = MLpars1[4], K_0 = MLpars1[5], z = MLpars1[6], gamma_0 = MLpars1[7], alpha = MLpars1[8], lambda_a0 = MLpars1[9], beta = MLpars1[10], loglik = ML, df = length(initparsopt), conv = unlist(out$conv))
    }
  }
  print_parameters_and_loglik(pars = MLpars1,
                              loglik = ML,
                              verbose = verbose,
                              type = 'global_ML',
                              distance_dep = distance_dep)
  return(invisible(out2))
}
