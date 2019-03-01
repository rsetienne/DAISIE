odeproc <- function(
  probs,
  times,
  fun,
  params,
  rtol,
  atol,
  method
)
{
  tshift = -abs(params[11])
  params1 = c(params[1:5],params[12:13])
  params2 = c(params[6:10],params[12:13])
  if(times[1] < tshift & times[2] < tshift)
  {
    #y = deSolve::ode(probs,times[1:2],fun,params1,rtol = rtol,atol = atol,method = method)
    y = DAISIE_integrate(probs,times[1:2],fun,params1,rtol = rtol,atol = atol,method = method)
  } else
    if(times[1] > tshift & times[2] > tshift)
    {
      #y = deSolve::ode(probs,times[1:2],fun,params2,rtol = rtol,atol = atol,method = metho d)
      y = DAISIE_integrate(probs,times[1:2],fun,params2,rtol = rtol,atol = atol,method = method)
    } else
      if(times[1] < tshift & times[2] > tshift)
      {
        #y = deSolve::ode(probs,c(times[1],tshift),fun,params1,rtol = rtol,atol = atol,method = method)
        y = DAISIE_integrate(probs,c(times[1],tshift),fun,params1,rtol = rtol,atol = atol,method = method)
        probs = y[2,2:ncol(y)]
        #y = deSolve::ode(probs,c(tshift,times[2]),fun,params2,rtol = rtol,atol = atol,method = method)
        y = DAISIE_integrate(probs,c(tshift,times[2]),fun,params2,rtol = rtol,atol = atol,method = method)
      }
  return(y)
}

divdepvecproc = function(
  params,
  lx,
  k1,
  ddep,
  times,
  type
)
{
  tshift = -abs(params[11])
  if(type == 'col')
  {
    a = 3
  } else
  {
    a = 0
  } 
  if(times < tshift)
  {
    return(divdepvec(params[1 + a],params[3],lx,k1,ddep))
  } else
  {
    return(divdepvec(params[6 + a],params[8],lx,k1,ddep))
  }
}       

DAISIE_SR_loglik_CS_M1 <- DAISIE_SR_loglik <- function(
  pars1,
  pars2,
  brts,             
  stac,
  missnumspec,
  methode = "lsodes",
  abstolint = abstolint,
  reltolint = reltolint
)
{
  # brts = branching times (positive, from present to past)
  # - max(brts) = age of the island
  # - next largest brts = stem age / time of divergence from the mainland
  # The interpretation of this depends on stac (see below)
  # For stac = 0, there is no other value.
  # For stac = 1 and stac = 5, this is the time since divergence from the immigrant's sister on the mainland.
  # The immigrant must have immigrated at some point since then.
  # For stac = 2 and stac = 3, this is the time since divergence from the mainland.
  # The immigrant that established the clade on the island must have immigrated precisely at this point.
  # For stac = 3, it must have reimmigrated, but only after the first immigrant had undergone speciation.
  # - min(brts) = most recent branching time (only for stac = 2, or stac = 3)
  # pars1 = model parameters
  # - pars1[1] = lac = (initial) cladogenesis rate
  # - pars1[2] = mu = extinction rate
  # - pars1[3] = K = maximum number of species possible in the clade
  # - pars1[4] = gam = (initial) immigration rate
  # - pars1[5] = laa = (initial) anagenesis rate
  # pars2 = model settings
  # - pars2[1] = lx = length of ODE variable x
  # - pars2[2] = ddep = diversity-dependent model,mode of diversity-dependence
  #  . ddep == 0 : no diversity-dependence
  #  . ddep == 1 : linear dependence in speciation rate (anagenesis and cladogenesis)
  #  . ddep == 11 : linear dependence in speciation rate and immigration rate
  #  . ddep == 3 : linear dependence in extinction rate
  # - pars2[3] = cond = conditioning
  #  . cond == 0 : no conditioning
  #  . cond == 1 : conditioning on presence on the island (not used in this single loglikelihood)
  # - pars2[4] = parameters and likelihood should be printed (1) or not (0)
  # stac = status of the clade formed by the immigrant
  #  . stac == 0 : immigrant is not present and has not formed an extant clade
  #  . stac == 1 : immigrant is present but has not formed an extant clade
  #  . stac == 2 : immigrant is not present but has formed an extant clade
  #  . stac == 3 : immigrant is present and has formed an extant clade
  #  . stac == 4 : immigrant is present but has not formed an extant clade, and it is known when it immigrated.
  #  . stac == 5 : immigrant is not present and has not formed an extant clade, but only an endemic species
  # missnumspec = number of missing species
  
  if(is.na(pars2[4]))
  {
    pars2[4] = 0
  }
  ddep = pars2[2]
  cond = pars2[3]
  if(cond > 0)
  {
    cat("Conditioning has not been implemented and may not make sense. Cond is set to 0.\n")
  }
  
  lac = pars1[1]
  mu = pars1[2]
  K = pars1[3]
  if(ddep == 0)
  {
    K = Inf
  }
  #gam = pars1[4]
  #laa = pars1[5]
  lac2 = pars1[6]
  mu2 = pars1[7]
  K2 = pars1[8]
  if(ddep == 0)
  {
    K2 = Inf
  }
  #gam2 = pars1[9]
  #laa2 = pars1[10]
  
  abstol = 1e-16
  reltol = 1e-10
  brts = -sort(abs(as.numeric(brts)),decreasing = TRUE)
  tshift = -abs(pars1[11])
  if(length(brts) == 1 & sum(brts == 0) == 1)
  {
    stop('The branching times contain only a 0. This means the island emerged at the present which is not allowed.');
    loglik = -Inf
    return(loglik)
  }
  if(sum(brts == 0) == 0)
  {
    brts[length(brts) + 1] = 0
  }
  # for stac = 0 and stac = 1, brts will contain origin of island and 0; length = 2; no. species should be 0
  # for stac = 1, brts will contain origin of island and 0; length = 2; no. species should be 1
  # for stac = 4, brts will contain origin of island, colonization event and 0; length = 3; no. species should be 1
  # for stac = 2, brts will contain origin of island, colonization event, branching times, 0; no. species should be no. branching times + 1
  # for stac = 3, brts will contain origin of island, colonization event, branching times, 0; no. species should be no. branching times + 2
  # for stac = 5, brts will contain origin of island, and 0; length = 2; number of species should be 1
  S = 0 * (stac == 0) + (stac == 1 || stac == 4 || stac == 5) + (length(brts) - 2) * (stac == 2) + (length(brts) - 1) * (stac == 3) + (length(brts) - 1) * (stac == 6) + length(brts) * (stac == 7)
  #S = length(brts) - (stac %% 2 == 1) - 2 * (stac %% 2 == 0) # old code before introduction of stac 6 and 7
  S2 = S - (stac == 1) - (stac == 3) - (stac == 4) - (stac == 7)
  loglik = -lgamma(S2 + missnumspec + 1) + lgamma(S2 + 1) + lgamma(missnumspec + 1)
  if(min(pars1) < 0)
  {
    cat('One or more parameters are negative.\n')
    loglik = -Inf
    return(loglik)
  }
  kshift = length(which(brts < tshift)) + 1 - (stac %% 2 == 1) - 2 * (stac %% 2 == 0) 
  if((ddep == 1 | ddep == 11) & (ceiling(K) < kshift | ceiling(K2) < (S + missnumspec)))
  {
    cat('The proposed values of K and/or K2 are incompatible with the number of species in the clade. Likelihood for this parameter set will be set to -Inf.\n')
    loglik = -Inf
    return(loglik)
  }
  if(lac == Inf & lac2 == Inf & mu != Inf & mu2 != Inf & missnumspec == 0)
  {
    loglik = DAISIE_loglik_high_lambda(pars1,-brts,stac)
  } else {
    if(ddep == 1 | ddep == 11)
    {
      lx = min(1 + max(missnumspec,ceiling(K),ceiling(K2)),DDD::roundn(pars2[1]) + missnumspec)
    } else {
      lx = DDD::roundn(pars2[1]) + missnumspec
    }
    if(loglik > -Inf)
    { 
      # in all cases we integrate from the origin of the island to the first branching point (stac > 1) or to the present (stac <= 1)
      probs = rep(0,2 * lx + 1)
      probs[1] = 1
      k1 = 0
      y = odeproc(probs,brts[1:2],DAISIE_loglik_rhs,c(pars1,k1,ddep),rtol = reltolint,atol = abstolint,method = methode)
      probs = y[2,2:(2 * lx + 2)]
      cp = checkprobs(lv = 2 * lx,loglik,probs); loglik = cp[[1]]; probs = cp[[2]]      
      if(stac == 0)
        # for stac = 0, the integration is from the origin of the island until the present
        # and we evaluate the probability of no clade being present and no immigrant species,
        # but there can be missing species
      {     
        loglik = loglik + log(probs[1 + missnumspec])
      } else {
        if(stac == 1 || stac == 5)
          # for stac = 1, the integration is from the maximum colonization time (usually the
          # island age + tiny time unit) until the present, where we set all probabilities where
          # the immigrant is already present to 0
          # and we evaluate the probability of the immigrant species being present,
          # but there can be missing species 
          # for stac = 5, we do exactly the same, but we evaluate the probability of an endemic species being present alone.          
        {         
          probs[(lx + 1):(2 * lx)] = 0
          y = odeproc(probs,brts[2:3],DAISIE_loglik_rhs,c(pars1,k1,ddep),rtol = reltolint,atol = abstolint,method = methode)
          probs = y[2,2:(2 * lx + 2)]
          cp = checkprobs(lv = 2 * lx,loglik,probs); loglik = cp[[1]]; probs = cp[[2]]               
          loglik = loglik + log(probs[(stac == 1) * lx + (stac == 5) + 1 + missnumspec])
        } else {
          # for stac > 1, but not 5, integration is then from the colonization event until the first branching time (stac = 2 and 3) or the present (stac = 4). We add a set of equations for Q_M,n, the probability that the process is compatible with the data, and speciation has not happened; during this time immigration is not allowed because it would alter the colonization time. After speciation, colonization is allowed again (re-immigration)
          # all probabilities of states with the immigrant present are set to zero and all probabilities of states with endemics present are transported to the state with the colonist present waiting for speciation to happen. We also multiply by the (possibly diversity-dependent) immigration rate
          if(stac == 6 || stac == 7)
          {
            probs[(lx + 1):(2 * lx)] = 0
            y = odeproc(probs,brts[2:3],DAISIE_loglik_rhs,c(pars1,k1,ddep),rtol = reltolint,atol = abstolint,method = methode)
            probs = y[2,2:(2 * lx + 2)]
            cp = checkprobs(lv = 2 * lx,loglik,probs); loglik = cp[[1]]; probs = cp[[2]] 
            k1 = 1
          }
          if(stac == 2 || stac == 3 || stac == 4)
          {
            gamvec = divdepvecproc(pars1,lx,k1,ddep * (ddep == 11 | ddep == 21),brts[2],'col')
            probs[(2 * lx + 1):(3 * lx)] = gamvec[1:lx] * probs[1:lx]
            probs[1:(2 * lx)] = 0        
            k1 = 1
            y = odeproc(probs,brts[2:3],DAISIE_loglik_rhs2,c(pars1,k1,ddep),rtol = reltolint,atol = abstolint,method = methode)
            probs = y[2,2:(3 * lx + 1)]
            cp = checkprobs2(lx,loglik,probs); loglik = cp[[1]]; probs = cp[[2]]
          }
          if(stac == 4)
            # if stac = 4, we're done and we take an element from Q_M,n
          {
            loglik = loglik + log(probs[2 * lx + 1 + missnumspec])
          } else {         
            # for stac = 2 and 3, at the first branching point all probabilities of states Q_M,n are transferred to probabilities where only endemics are present. Then go through the branching points.
            S1 = length(brts) - 1
            startk = (3 - (stac == 6 || stac == 7))
            if(S1 >= startk)
            {
              lacvec = divdepvecproc(pars1,lx,k1,ddep,brts[3],'spec')
              if(stac == 2 || stac == 3)
              {
                probs[1:lx] = lacvec[1:lx] * (probs[1:lx] + probs[(2 * lx + 1):(3 * lx)])
                probs[(lx + 1):(2 * lx)] = lacvec[2:(lx + 1)] * probs[(lx + 1):(2 * lx)]
                probs = probs[-c((2 * lx + 2):(3 * lx))]
                probs[2 * lx + 1] = 0
              }
              if(stac == 6 || stac == 7)
              {
                probs2 = rep(0,2 * lx + 1)
                probs2[(1:(lx - 1))] = probs[(2:lx)] + 1/(2:lx) * probs[(lx + 1):(2 * lx - 1)]
                probs2[lx] = 1/(lx + 1) * probs[2 * lx]
                probs2[(lx + 1):(2 * lx - 1)] = (1:(lx - 1))/(2:lx) * probs[(lx + 2):(2 * lx)]
                probs = probs2             
                rm(probs2)  
                probs[1:lx] = lacvec[1:lx] * probs[1:lx]
                probs[(lx + 1):(2 * lx)] = lacvec[2:(lx + 1)] * probs[(lx + 1):(2 * lx)]
              }
              for(k in startk:S1)
              {
                k1 = k - (stac == 2 || stac == 3)
                y = odeproc(probs,brts[k:(k+1)],DAISIE_loglik_rhs,c(pars1,k1,ddep),rtol = reltolint,atol = abstolint,method = methode)
                probs = y[2,2:(2 * lx + 2)]
                cp = checkprobs2(lx,loglik,probs); loglik = cp[[1]]; probs = cp[[2]]
                if(k < S1)
                {
                  # speciation event      
                  lacvec = divdepvecproc(pars1,lx,k1,ddep,brts[k + 1],'spec')
                  probs[1:(2 * lx)] = c(lacvec[1:lx],lacvec[2:(lx + 1)]) * probs[1:(2 * lx)]
                }
              }            
            }
            # we evaluate the probability of the phylogeny with any missing species at the present without (stac = 2) or with (stac = 3) the immigrant species; there can be no missing species for stac = 4
            loglik = loglik + log(probs[(stac == 3 || stac == 7) * lx + 1 + missnumspec])
          }   
        }     
      }           
    }
  }
  #print(head(probs,n = 15))
  
  if(pars2[4] >= 1)
  {
    # if (length(pars1 > 5)) {
      s1 = sprintf('Status of colonist: %d, Parameters: %f %f %f %f %f ', stac, pars1[1], pars1[2], pars1[3], pars1[4], pars1[5])
    # }
    s2 = sprintf(', Loglikelihood: %f',loglik)
    cat(s1,s2,"\n",sep = "")
    utils::flush.console()
  }
  
  return(as.numeric(loglik))
}



#' Computes the loglikelihood of the DAISIE model with clade-specific
#' diversity-dependence given data and a set of model parameters that may shift
#' at some time
#' 
#' Computes the loglikelihood of the DAISIE model with clade-specific
#' diversity-dependence given colonization and branching times for lineages on
#' an island, and a set of model parameters that may shift at some time
#' 
#' The output is a loglikelihood value
#' 
#' @aliases DAISIE_SR_loglik_CS DAISIE_SR_loglik_all
#' @param pars1 Contains the model parameters: \cr \cr \code{pars1[1]}
#' corresponds to lambda^c (cladogenesis rate) \cr \code{pars1[2]} corresponds
#' to mu (extinction rate) \cr \code{pars1[3]} corresponds to K (clade-level
#' carrying capacity) \cr \code{pars1[4]} corresponds to gamma (immigration
#' rate) \cr \code{pars1[5]} corresponds to lambda^a (anagenesis rate) \cr
#' \code{pars1[6]} corresponds to lambda^c (cladogenesis rate) after the shift
#' \cr \code{pars1[7]} corresponds to mu (extinction rate) after the shift \cr
#' \code{pars1[8]} corresponds to K (clade-level carrying capacity) after the
#' shift \cr \code{pars1[9]} corresponds to gamma (immigration rate) after the
#' shift \cr \code{pars1[10]} corresponds to lambda^a (anagenesis rate) after
#' the shift \cr \code{pars1[11]} corresponds to the time of shift \cr
#' @param pars2 Contains the model settings \cr \cr \code{pars2[1]} corresponds
#' to lx = length of ODE variable x \cr \code{pars2[2]} corresponds to ddmodel
#' = diversity-dependent model, model of diversity-dependence, which can be one
#' of\cr \cr ddmodel = 0 : no diversity dependence \cr ddmodel = 1 : linear
#' dependence in speciation rate \cr ddmodel = 11: linear dependence in
#' speciation rate and in immigration rate \cr ddmodel = 2 : exponential
#' dependence in speciation rate\cr ddmodel = 21: exponential dependence in
#' speciation rate and in immigration rate\cr \cr \code{pars2[3]} corresponds
#' to cond = setting of conditioning\cr \cr cond = 0 : conditioning on island
#' age \cr cond = 1 : conditioning on island age and non-extinction of the
#' island biota \cr \cr \code{pars2[4]} sets whether parameters and likelihood
#' should be printed (1) or not (0)
#' @param datalist Data object containing information on colonisation and
#' branching times. This object can be generated using the DAISIE_dataprep
#' function, which converts a user-specified data table into a data object, but
#' the object can of course also be entered directly. It is an R list object
#' with the following elements.\cr The first element of the list has two or
#' three components: \cr \cr \code{$island_age} - the island age \cr Then,
#' depending on whether a distinction between types is made, we have:\cr
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
#' @param methode Method of the ODE-solver. See package deSolve for details.
#' Default is "lsodes"
#' @param CS_version For internal testing purposes only. Default is 1, the
#' original DAISIE code.
#' @param abstolint Absolute tolerance of the integration
#' @param reltolint Relative tolerance of the integration
#' @return The loglikelihood
#' @author Rampal S. Etienne & Bart Haegeman
#' @seealso \code{\link{DAISIE_ML}}, \code{\link{DAISIE_sim}}
#' @references Valente, L.M., A.B. Phillimore and R.S. Etienne (2015).
#' Equilibrium and non-equilibrium dynamics simultaneously operate in the
#' Galapagos islands. Ecology Letters 18: 844-852.
#' @keywords models
#' @examples
#' 
#' utils::data(Galapagos_datalist_2types)
#' pars1 = c(0.195442017,0.087959583,Inf,0.002247364,0.873605049,
#'           3755.202241,8.909285094,14.99999923,0.002247364,0.873605049,0.163)
#' pars2 = c(100,11,0,1)
#' DAISIE_loglik_all(pars1,pars2,Galapagos_datalist_2types)
#' 
#' @export DAISIE_SR_loglik_CS
#' @export DAISIE_SR_loglik_all
DAISIE_SR_loglik_CS <- DAISIE_SR_loglik_all <- function(
  pars1,
  pars2,
  datalist,
  methode = "lsodes",
  CS_version = 1,
  abstolint = 1E-16,
  reltolint = 1E-10
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
  #  . stac == 5 : immigrant is not present and has not formed an extant clade, but only an endemic species
  # datalist[[,]][3] = list with number of missing species in clades for stac = 2 and stac = 3;
  # for stac = 0 and stac = 1, this number equals 0.
  # pars1 = model parameters
  # - pars1[1] = lac = (initial) cladogenesis rate
  # - pars1[2] = mu = extinction rate
  # - pars1[3] = K = maximum number of species possible in the clade
  # - pars1[4] = gam = (initial) immigration rate
  # - pars1[5] = laa = (initial) anagenesis rate
  # - pars1[6]...pars1[10] = same as pars1[1]...pars1[5], but after the shift
  # - pars1[11] = time of shift
  # pars2 = model settings
  # - pars2[1] = lx = length of ODE variable x
  # - pars2[2] = ddep = diversity-dependent model,mode of diversity-dependence
  #  . ddep == 0 : no diversity-dependence
  #  . ddep == 1 : linear dependence in speciation rate (anagenesis and cladogenesis)
  #  . ddep == 11 : linear dependence in speciation rate and immigration rate
  #  . ddep == 3 : linear dependence in extinction rate
  # - pars2[3] = cond = conditioning
  #  . cond == 0 : no conditioning
  #  . cond == 1 : conditioning on presence on the island (not used in this single loglikelihood)
  # - pars2[4] = parameters and likelihood should be printed (1) or not (0)
  
  pars1 = as.numeric(pars1)
  cond = pars2[3]
  logp0 = DAISIE_SR_loglik_CS_M1(pars1 = pars1,pars2 = pars2,brts = datalist[[1]]$island_age,stac = 0,missnumspec = 0,methode = methode,abstolint = abstolint,reltolint = reltolint)
  if(is.null(datalist[[1]]$not_present))
  {
    not_present = (datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2)
  } else {
    not_present = datalist[[1]]$not_present
  }
  loglik = not_present * logp0
  numimm = not_present + length(datalist) - 1
  logcond = (cond == 1) * log(1 - exp(numimm * logp0))
  loglik = loglik - logcond
  if(length(datalist) > 1)
  {
    for(i in 2:length(datalist))
    {
      loglik = loglik + DAISIE_SR_loglik_CS_M1(pars1 = pars1,pars2 = pars2,brts = datalist[[i]]$branching_times,stac = datalist[[i]]$stac,missnumspec = datalist[[i]]$missing_species,methode = methode,abstolint = abstolint,reltolint = reltolint)
    }
  }
  return(loglik)
}
