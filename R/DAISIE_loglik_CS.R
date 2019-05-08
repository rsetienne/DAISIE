DAISIE_loglik_rhs_precomp <- function(pars,lx)
{
  lac = pars[1]
  mu = pars[2]
  K = pars[3]
  gam = pars[4]
  laa = pars[5]
  kk = pars[6]
  ddep = pars[7]
  
  
  nn = -2:(lx+2*kk+1)
  lnn = length(nn)
  nn = pmax(rep(0,lnn),nn)
  
  if(ddep == 0)
  {
    laavec = laa * rep(1,lnn)
    lacvec = lac * rep(1,lnn)
    muvec = mu * rep(1,lnn)
    gamvec = gam * rep(1,lnn)
  } else if(ddep == 1)
  {
    laavec = laa * rep(1,lnn)
    lacvec = pmax(rep(0,lnn),lac * (1 - nn/K))
    muvec = mu * rep(1,lnn)
    gamvec = gam * rep(1,lnn)
  } else if(ddep == 2)
  {
    laavec = laa * rep(1,lnn)
    lacvec = pmax(rep(0,lnn),lac * exp(-nn/K))
    muvec = mu * rep(1,lnn)
    gamvec = gam * rep(1,lnn)
  } else if(ddep == 11)
  {
    laavec = laa * rep(1,lnn)
    lacvec = pmax(rep(0,lnn),lac * (1 - nn/K))
    muvec = mu * rep(1,lnn)
    gamvec = pmax(rep(0,lnn),gam * (1 - nn/K))
  } else if(ddep == 21)
  {
    laavec = laa * rep(1,lnn)
    lacvec = pmax(rep(0,lnn),lac * exp(-nn/K))
    muvec = mu * rep(1,lnn)
    gamvec = pmax(rep(0,lnn),gam * exp(-nn/K))
  } else if(ddep == 3)
  {
    laavec = laa * rep(1,lnn)
    lacvec = lac * rep(1,lnn)
    muvec = mu * (1 + nn/K)
    gamvec = gam * rep(1,lnn)
  }
  return(c(laavec,lacvec,muvec,gamvec,nn,kk))
}

DAISIE_loglik_rhs = function(t,x,parsvec)
{
  kk <- parsvec[length(parsvec)]
  lx <- (length(x) - 1)/2
  lnn <- lx + 4 + 2 * kk
  laavec <- parsvec[1:lnn]
  lacvec <- parsvec[(lnn + 1):(2 * lnn)]
  muvec <- parsvec[(2 * lnn + 1):(3 * lnn)]
  gamvec <- parsvec[(3 * lnn + 1):(4 * lnn)]
  nn <- parsvec[(4 * lnn + 1):(5 * lnn)]
  
  xx1 = c(0,0,x[1:lx],0)
  xx2 = c(0,0,x[(lx + 1):(2 * lx)],0)
  xx3 = x[2 * lx + 1]
  
  nil2lx = 3:(lx + 2)
  
  il1 = nil2lx+kk-1
  il2 = nil2lx+kk+1
  il3 = nil2lx+kk
  il4 = nil2lx+kk-2
  
  in1 = nil2lx+2*kk-1
  in2 = nil2lx+1
  in3 = nil2lx+kk
  
  ix1 = nil2lx-1
  ix2 = nil2lx+1
  ix3 = nil2lx
  ix4 = nil2lx-2
  
  dx1 = laavec[il1 + 1] * xx2[ix1] + lacvec[il4 + 1] * xx2[ix4] + muvec[il2 + 1] * xx2[ix3] +
    lacvec[il1] * nn[in1] * xx1[ix1] + muvec[il2] * nn[in2] * xx1[ix2] +
    -(muvec[il3] + lacvec[il3]) * nn[in3] * xx1[ix3] +
    -gamvec[il3] * xx1[ix3]
  dx1[1] = dx1[1] + laavec[il3[1]] * xx3 * (kk == 1)
  dx1[2] = dx1[2] + 2 * lacvec[il3[1]] * xx3 * (kk == 1)
  
  dx2 = gamvec[il3] * xx1[ix3] +
    lacvec[il1 + 1] * nn[in1] * xx2[ix1] + muvec[il2 + 1] * nn[in2] * xx2[ix2] +
    -(muvec[il3 + 1] + lacvec[il3 + 1]) * nn[in3 + 1] * xx2[ix3] +
    -laavec[il3 + 1] * xx2[ix3]
  
  dx3 = -(laavec[il3[1]] + lacvec[il3[1]] + gamvec[il3[1]] + muvec[il3[1]]) * xx3
  
  return(list(c(dx1,dx2,dx3)))
}

DAISIE_loglik_rhs2 = function(t,x,parsvec)
{
  kk <- parsvec[length(parsvec)]
  lx <- (length(x))/3
  lnn <- lx + 4 + 2 * kk
  laavec <- parsvec[1:lnn]
  lacvec <- parsvec[(lnn + 1):(2 * lnn)]
  muvec <- parsvec[(2 * lnn + 1):(3 * lnn)]
  gamvec <- parsvec[(3 * lnn + 1):(4 * lnn)]
  nn <- parsvec[(4 * lnn + 1):(5 * lnn)]
  
  xx1 = c(0,0,x[1:lx],0)
  xx2 = c(0,0,x[(lx + 1):(2 * lx)],0)
  xx3 = c(0,0,x[(2 * lx + 1):(3 * lx)],0)
  
  nil2lx = 3:(lx + 2)
  
  il1 = nil2lx+kk-1
  il2 = nil2lx+kk+1
  il3 = nil2lx+kk
  il4 = nil2lx+kk-2
  
  in1 = nil2lx+2*kk-1
  in2 = nil2lx+1
  in3 = nil2lx+kk
  in4 = nil2lx-1
  
  ix1 = nil2lx-1
  ix2 = nil2lx+1
  ix3 = nil2lx
  ix4 = nil2lx-2
  
  # inflow:
  # anagenesis in colonist when k = 1: Q_M,n -> Q^1_n; n+k species present
  # cladogenesis in colonist when k = 1: Q_M,n-1 -> Q^1_n; n+k-1 species present; rate twice
  # anagenesis of reimmigrant: Q^M,k_n-1 -> Q^k,n; n+k-1+1 species present
  # cladogenesis of reimmigrant: Q^M,k_n-2 -> Q^k,n; n+k-2+1 species present; rate once
  # extinction of reimmigrant: Q^M,k_n -> Q^k,n; n+k+1 species present
  # cladogenesis in one of the n+k-1 species: Q^k_n-1 -> Q^k_n; n+k-1 species present; rate twice for k species
  # extinction in one of the n+1 species: Q^k_n+1 -> Q^k_n; n+k+1 species present
  # outflow:
  # all events with n+k species present
  dx1 = (laavec[il3] * xx3[ix3] + 2 * lacvec[il1] * xx3[ix1]) * (kk == 1) + 
    laavec[il1 + 1] * xx2[ix1] + lacvec[il4 + 1] * xx2[ix4] + muvec[il2 + 1] * xx2[ix3] +
    lacvec[il1] * nn[in1] * xx1[ix1] + muvec[il2] * nn[in2] * xx1[ix2] +
    -(muvec[il3] + lacvec[il3]) * nn[in3] * xx1[ix3] - gamvec[il3] * xx1[ix3]
  
  # inflow:
  # immigration when there are n+k species: Q^k,n -> Q^M,k_n; n+k species present
  # cladogenesis in n+k-1 species: Q^M,k_n-1 -> Q^M,k_n; n+k-1+1 species present; rate twice for k species
  # extinction in n+1 species: Q^M,k_n+1 -> Q^M,k_n; n+k+1+1 species present
  # outflow:
  # all events with n+k+1 species present
  dx2 = gamvec[il3] * xx1[ix3] +
    lacvec[il1 + 1] * nn[in1] * xx2[ix1] + muvec[il2 + 1] * nn[in2] * xx2[ix2] +
    -(muvec[il3 + 1] + lacvec[il3 + 1]) * nn[in3 + 1] * xx2[ix3] +
    -laavec[il3 + 1] * xx2[ix3]
  
  # only when k = 1         
  # inflow:
  # cladogenesis in one of the n-1 species: Q_M,n-1 -> Q_M,n; n+k-1 species present; rate once
  # extinction in one of the n+1 species: Q_M,n+1 -> Q_M,n; n+k+1 species present
  # outflow:
  # all events with n+k species present
  dx3 = lacvec[il1] * nn[in4] * xx3[ix1] + muvec[il2] * nn[in2] * xx3[ix2] +
    -(lacvec[il3] + muvec[il3]) * nn[in3] * xx3[ix3] +
    -(laavec[il3] + gamvec[il3]) * xx3[ix3]
  
  return(list(c(dx1,dx2,dx3)))
}

checkprobs = function(lv,loglik,probs)
{
  probs = probs * (probs > 0)
  if(is.na(sum(probs[1:lv])) || is.nan(sum(probs)))
  {
    cat('Numerical issues encountered\n')
    loglik = -Inf
  } else if(sum(probs[1:lv]) <= 0)
  {
    cat('Numerical issues encountered\n')
    loglik = -Inf
  } else {
    loglik = loglik + log(sum(probs[1:lv]))
    probs[1:lv] = probs[1:lv]/sum(probs[1:lv])
  }
  return(list(loglik,probs))
}

checkprobs2 = function(lx,loglik,probs)
{
  probs = probs * (probs > 0)
  if(is.na(sum(probs)) || is.nan(sum(probs)))
  {
    cat('Numerical issues encountered\n')
    loglik = -Inf
  } else if(sum(probs) <= 0)
  {
    cat('Numerical issues encountered\n')
    loglik = -Inf
  } else {
    loglik = loglik + log(sum(probs))
    probs = probs/sum(probs)
  }   
  return(list(loglik,probs))
}

divdepvec <- function(lacgam,pars1,lx,k1,ddep,island_ontogeny = NA)
{
  if(!is.na(island_ontogeny))
  {
    lacgamK <- divdepvec_time(lacgam,pars1,lx,k1,ddep,island_ontogeny)
    lacgam <- lacgamK[1]
    K <- lacgamK[2]
  } else
  {
    K <- pars1[1]
  }
  return(divdepvec1(lacgam,K,lx,k1,ddep))
}

divdepvec1 = function(lacgam,K,lx,k1,ddep)
{
  if(ddep == 1 | ddep == 11)
  {
    vec = pmax(rep(0,lx + 1),lacgam * (1 - ((0:lx) + k1) / K))
  } else {
    if(ddep == 2 | ddep == 21)
    {
      vec = pmax(rep(0,lx + 1),lacgam * exp(-((0:lx)+k1) / K))
    } else {
      if(ddep == 0 | ddep == 3)
      {
        vec = lacgam * rep(1,lx + 1)
      }
    }
  }
  return(vec)
}

DAISIE_loglik_CS_M1 <- DAISIE_loglik <- function(
  pars1,
  pars2,
  brts,
  stac,
  missnumspec,
  methode = "lsodes",
  abstolint = 1E-16,
  reltolint = 1E-10
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
  # - pars2[5] = island ontonogeny. If NULL, then constant ontogeny is assumed
  # missnumspec = number of missing species
  # stac = status of the clade formed by the immigrant
  #  . stac == 1 : immigrant is present but has not formed an extant clade
  #  . stac == 2 : immigrant is not present but has formed an extant clade
  #  . stac == 3 : immigrant is present and has formed an extant clade
  #  . stac == 4 : immigrant is present but has not formed an extant clade, and it is known when it immigrated.
  #  . stac == 5 : immigrant is not present and has not formed an extant clade, but only an endemic species
  #  . stac == 6 : like 2, but with max colonization time
  #  . stac == 7 : like 3, but with max colonization time

  # Stop laa from being inf and return -Inf  
  if (is.infinite(pars1[5])) {
    return(-Inf)
  }
  
  if(is.na(pars2[4]))
  {
    pars2[4] = 0
  }
  ddep = pars2[2]
  cond = pars2[3]
  # TODO: check if pars2[5] should be NA of if this never happens
  # if (is.na(pars2[5])) { 
  #   pars2[5] <- 0
  # }
  island_ontogeny <- pars2[5]
  if(cond > 0)
  {
    cat("Conditioning has not been implemented and may not make sense. Cond is set to 0.\n")
  }
  
  if (is.na(island_ontogeny)) # This calls the old code that doesn't expect 
    # ontogeny
  {
    lac = pars1[1]
    mu = pars1[2]
    K = pars1[3]
    if(ddep == 0)
    {
      K = Inf
    }
    gam = pars1[4]
    laa = pars1[5]
    pars1_in_divdepvec_call <- K
  } else {
    #pars1[1:4] = Apars
    #pars1[5] = lac0
    #pars1[6:7] = mupars
    #pars1[8] = K0
    #pars1[9] = gam0
    #pars1[10] = laa
    #pars1[11] = island_ontogeny
    pars1[11] <- island_ontogeny
    
    if (pars1[11] == 0 && pars1[6] != pars1[7]) {
      warning("mu_min and mu_max are not equal! Setting mu_max = mu_min")
      pars1[7] <- pars1[6]
    }
    
    lac <- as.numeric(pars1[5])
    K <- as.numeric(pars1[8])
    gam <- as.numeric(pars1[9])
    pars1_in_divdepvec_call <- pars1
  }
  
  brts = -sort(abs(as.numeric(brts)),decreasing = TRUE)
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
  # for stac = 0, brts will contain origin of island and 0; length = 2; no. species should be 0
  # for stac = 1, brts will contain origin of island, maximum colonization time (usually island age) and 0; length = 3; no. species should be 1
  # for stac = 2, brts will contain origin of island, colonization event, branching times, 0; no. species should be no. branching times + 1
  # for stac = 3, brts will contain origin of island, colonization event, branching times, 0; no. species should be no. branching times + 2
  # for stac = 4, brts will contain origin of island, colonization event and 0; length = 3; no. species should be 1
  # for stac = 5, brts will contain origin of island, maximum colonization time (usually island age), and 0; length = 2; number of species should be 1
  # for stac = 6, brts will contain origin of island, maximum colonization time (usually island age), branching times and 0; number of species should be no. branching times + 1
  # for stac = 7, brts will contain origin of island, maximum colonization time (usually island age), branching times and 0; number of species should be no. branching times + 2
  S = 0 * (stac == 0) + (stac == 1 || stac == 4 || stac == 5) + (length(brts) - 2) * (stac == 2) + (length(brts) - 1) * (stac == 3) + (length(brts) - 2) * (stac == 6) + (length(brts) - 1) * (stac == 7)
  #S = length(brts) - (stac %% 2 == 1) - 2 * (stac %% 2 == 0) # old code before introduction of stac 6 and 7
  S2 = S - (stac == 1) - (stac == 3) - (stac == 4) - (stac == 7)
  loglik = -lgamma(S2 + missnumspec + 1) + lgamma(S2 + 1) + lgamma(missnumspec + 1)
  if(min(pars1) < 0)
  {
    cat('One or more parameters are negative.\n')
    loglik = -Inf
    return(loglik)
  }
  if((ddep == 1 | ddep == 11) & ceiling(K) < (S + missnumspec))
  {
    cat('The proposed value of K is incompatible with the number of species in the clade. Likelihood for this parameter set will be set to -Inf.\n')
    loglik = -Inf
    return(loglik)
  }
  if(lac == Inf & missnumspec == 0 & length(pars1) == 5)
  {
    loglik = DAISIE_loglik_high_lambda(pars1,-brts,stac)
  } else {
    if(ddep == 1 | ddep == 11)
    {
      lx = min(1 + max(missnumspec,ceiling(K)),DDD::roundn(pars2[1]) + missnumspec)
    } else {
      lx = DDD::roundn(pars2[1]) + missnumspec
    }
    if(loglik > -Inf)
    { 
      # in all cases we integrate from the origin of the island to the colonization event (stac 2, 3, 4), the first branching point (stac = 6, 7) or to the present (stac = 0, 1, 5)
      probs = rep(0,2 * lx + 1)
      probs[1] = 1
      k1 = 0
      #y = deSolve::ode(probs,brts[1:2],DAISIE_loglik_rhs,c(pars1,k1,ddep),rtol = reltolint,atol = abstolint,method = methode)
      y = DAISIE_integrate(probs,brts[1:2],DAISIE_loglik_rhs,c(pars1,k1,ddep),rtol = reltolint,atol = abstolint,method = methode)
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
          #y = deSolve::ode(probs,brts[2:3],DAISIE_loglik_rhs,c(pars1,k1,ddep),rtol = reltolint,atol = abstolint,method = methode)
          y = DAISIE_integrate(probs,brts[2:3],DAISIE_loglik_rhs,c(pars1,k1,ddep),rtol = reltolint,atol = abstolint,method = methode)
          probs = y[2,2:(2 * lx + 2)]
          cp = checkprobs(lv = 2 * lx,loglik,probs); loglik = cp[[1]]; probs = cp[[2]]               
          loglik = loglik + log(probs[(stac == 1) * lx + (stac == 5) + 1 + missnumspec])
        } else {
          # for stac = 2, 3, 4, integration is then from the colonization event until the first branching time (stac = 2 and 3) or the present (stac = 4). We add a set of equations for Q_M,n, the probability that the process is compatible with the data, and speciation has not happened; during this time immigration is not allowed because it would alter the colonization time. After speciation, colonization is allowed again (re-immigration)
          # all probabilities of states with the immigrant present are set to zero and all probabilities of states with endemics present are transported to the state with the colonist present waiting for speciation to happen. We also multiply by the (possibly diversity-dependent) immigration rate
          # for stac = 6 and 7, integration is from the maximum colonization time until the first branching time
          if(stac == 6 || stac == 7)
          {
            probs[(lx + 1):(2 * lx)] = 0
            #y = deSolve::ode(probs,brts[2:3],DAISIE_loglik_rhs,c(pars1,k1,ddep),rtol = reltolint,atol = abstolint,method = methode)
            y = DAISIE_integrate(probs,brts[2:3],DAISIE_loglik_rhs,c(pars1,k1,ddep),rtol = reltolint,atol = abstolint,method = methode)
            probs = y[2,2:(2 * lx + 2)]
            cp = checkprobs(lv = 2 * lx,loglik,probs); loglik = cp[[1]]; probs = cp[[2]] 
            k1 = 1
          }
          if(stac == 2 || stac == 3 || stac == 4)
          {
            t <- brts[2]
            gamvec = divdepvec(gam,c(pars1_in_divdepvec_call,t,0),lx,k1,ddep * (ddep == 11 | ddep == 21),island_ontogeny) # Problem may be here 30/3
            probs[(2 * lx + 1):(3 * lx)] = gamvec[1:lx] * probs[1:lx]
            probs[1:(2 * lx)] = 0        
            k1 = 1
            #y = deSolve::ode(probs,c(brts[2:3]),DAISIE_loglik_rhs2,c(pars1,k1,ddep),rtol = reltolint,atol = abstolint,method = methode)
            y = DAISIE_integrate(probs,c(brts[2:3]),DAISIE_loglik_rhs2,c(pars1,k1,ddep),rtol = reltolint,atol = abstolint,method = methode)
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
            startk = 3
            if(S1 >= startk)
            {
              t <- brts[startk]
              lacvec = divdepvec(lac,c(pars1_in_divdepvec_call,t,1),lx,k1,ddep,island_ontogeny)
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
                k1 = k - 1
                #y = deSolve::ode(probs,brts[k:(k+1)],DAISIE_loglik_rhs,c(pars1,k1,ddep),rtol = reltolint,atol = abstolint,method = methode)
                y = DAISIE_integrate(probs,brts[k:(k+1)],DAISIE_loglik_rhs,c(pars1,k1,ddep),rtol = reltolint,atol = abstolint,method = methode)
                probs = y[2,2:(2 * lx + 2)]
                cp = checkprobs2(lx,loglik,probs); loglik = cp[[1]]; probs = cp[[2]]
                if(k < S1)
                {
                  # speciation event      
                  t <- brts[k + 1]
                  lacvec = divdepvec(lac,c(pars1_in_divdepvec_call,t,1),lx,k1,ddep,island_ontogeny)
                  probs[1:(2 * lx)] = c(lacvec[1:lx],lacvec[2:(lx + 1)]) * probs[1:(2 * lx)]
                }
              }            
            }
            # we evaluate the probability of the phylogeny with any missing species at the present without (stac = 2 or stac = 6) or with (stac = 3 or stac = 7) the immigrant species
            loglik = loglik + log(probs[(stac == 3 || stac == 7) * lx + 1 + missnumspec])
          }   
        }     
      }           
    }
  }
  # print(head(probs,n = 15))
  
  if(pars2[4] >= 1)
  {
    if (length(pars1) == 11) { # CHANGE
      s1 = sprintf('Status of colonist: %d, Parameters: %f %f %f %f %f %f', stac, pars1[5], pars1[6], pars1[7], pars1[8], pars1[9], pars1[10])
    } else {
      s1 = sprintf('Status of colonist: %d, Parameters: %f %f %f %f %f ', stac, pars1[1], pars1[2], pars1[3], pars1[4], pars1[5])
    }
    s2 = sprintf(', Loglikelihood: %f',loglik)
    cat(s1,s2,"\n",sep = "")
    utils::flush.console()
  }
  
  return(as.numeric(loglik))
}

DAISIE_loglik_CS_choice = function(
  pars1,
  pars2,
  brts,
  stac,
  missnumspec,
  methode = "lsodes",
  CS_version = 1,
  abstolint = 1E-16,
  reltolint = 1E-10
)
{
  if(CS_version == 1)
  {
    loglik = DAISIE_loglik(pars1 = pars1,pars2 = pars2,brts = brts,stac = stac,missnumspec = missnumspec,methode = methode, abstolint = abstolint, reltolint = reltolint)
  } else
  {
    loglik = DAISIE_loglik_IW_M1(pars1 = pars1,pars2 = pars2,brts = brts,stac = stac,missnumspec = missnumspec,methode = methode, abstolint = abstolint, reltolint = reltolint)
  }
  return(loglik)
}


#' @name DAISIE_loglik_CS
#' @aliases DAISIE_loglik_all DAISIE_loglik_CS
#' @title Computes the loglikelihood of the DAISIE model with clade-specific
#' diversity-dependence given data and a set of model parameters
#' @description Computes the loglikelihood of the DAISIE model with clade-specific
#' diversity-dependence given colonization and branching times for lineages on
#' an island, and a set of model parameters. The output is a loglikelihood value
#' @param pars1 Contains the model parameters: \cr \cr \code{pars1[1]}
#' corresponds to lambda^c (cladogenesis rate) \cr \code{pars1[2]} corresponds
#' to mu (extinction rate) \cr \code{pars1[3]} corresponds to K (clade-level
#' carrying capacity) \cr \code{pars1[4]} corresponds to gamma (immigration
#' rate) \cr \code{pars1[5]} corresponds to lambda^a (anagenesis rate) \cr
#' \code{pars1[6]} corresponds to lambda^c (cladogenesis rate) for an optional
#' subset of the species \cr \code{pars1[7]} corresponds to mu (extinction
#' rate) for an optional subset of the species\cr \code{pars1[8]} corresponds
#' to K (clade-level carrying capacity) for an optional subset of the
#' species\cr \code{pars1[9]} corresponds to gamma (immigration rate) for an
#' optional subset of the species\cr \code{pars1[10]} corresponds to lambda^a
#' (anagenesis rate) for an optional subset of the species\cr \code{pars1[11]}
#' corresponds to p_f (fraction of mainland species that belongs to the second
#' subset of species\cr The elements 6:10 and 11 are optional, that is,pars1
#' should either contain 5, 10 or 11 elements. If 10, then the fraction of
#' potential colonists of type 2 is computed from the data. If 11, then
#' pars1[11] is used, overruling any information in the data.
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
#' on the island \cr or:\cr \code{$not_present_type1} - the number of mainland
#' lineages of type 1 that are not present on the island \cr
#' \code{$not_present_type2} - the number of mainland lineages of type 2 that
#' are not present on the island \cr \cr The remaining elements of the list
#' each contains information on a single colonist lineage on the island and has
#' 5 components:\cr \cr \code{$colonist_name} - the name of the species or
#' clade that colonized the island \cr \code{$branching_times} - island age and
#' stem age of the population/species in the case of Non-endemic,
#' Non-endemic_MaxAge and Endemic anagenetic species. For cladogenetic species
#' these should be island age and branching times of the radiation including
#' the stem age of the radiation.\cr \code{$stac} - the status of the colonist
#' \cr \cr * Non_endemic_MaxAge: 1 \cr * Endemic: 2 \cr * Endemic&Non_Endemic:
#' 3 \cr * Non_Endemic: 4 \cr * Endemic_Singleton_MaxAge: 5 \cr *
#' Endemic_Clade_MaxAge: 6 \cr * Endemic&Non_Endemic_Clade_MaxAge: 7 \cr \cr
#' \code{$missing_species} - number of island species that were not sampled for
#' particular clade (only applicable for endemic clades) \cr \code{$type1or2} -
#' whether the colonist belongs to type 1 or type 2 \cr
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
#' @export DAISIE_loglik_CS
#' @export DAISIE_loglik_all
DAISIE_loglik_CS <- DAISIE_loglik_all <- function(
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
  #  . stac == 6 : like 2, but with max colonization time
  #  . stac == 7 : like 3, but with max colonization time
  # datalist[[,]][3] = list with number of missing species in clades for stac = 2 and stac = 3;
  # for stac = 0 and stac = 1, this number equals 0.
  # pars1 = model parameters
  # - pars1[1] = lac = (initial) cladogenesis rate
  # - pars1[2] = mu = extinction rate
  # - pars1[3] = K = maximum number of species possible in the clade
  # - pars1[4] = gam = (initial) immigration rate
  # - pars1[5] = laa = (initial) anagenesis rate
  # - pars1[6]...pars1[10] = same as pars1[1]...pars1[5], but for a second type of immigrant
  # - pars1[11] = proportion of type 2 immigrants in species pool
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
  # - pars2[5] = island ontonogeny. If NA, then constant ontogeny is assumed
  
  pars1 = as.numeric(pars1)
  cond = pars2[3]
  endpars1 <- 5
  
  if(length(pars1) == 5 | !is.na(pars2[5])) # Normal no ont case
  {
    if(!is.na(pars2[5]))
    {
      endpars1 <- length(pars1)    
    }
    logp0 = DAISIE_loglik_CS_choice(pars1 = pars1,pars2 = pars2,brts = datalist[[1]]$island_age,stac = 0,missnumspec = 0,methode = methode,CS_version = CS_version,abstolint = abstolint,reltolint = reltolint)
    if(is.null(datalist[[1]]$not_present))
    {
      loglik = (datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2) * logp0
      numimm = (datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2) + length(datalist) - 1
    } else {
      loglik = datalist[[1]]$not_present * logp0
      numimm = datalist[[1]]$not_present + length(datalist) - 1
    }
    logcond = (cond == 1) * log(1 - exp(numimm * logp0))
    if(length(datalist) > 1)
    {
      for(i in 2:length(datalist))
      {
        datalist[[i]]$type1or2 = 1
      }
    }
  } else {
    numimm = datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2 + length(datalist) - 1
    numimm_type2 = length(which(unlist(datalist)[which(names(unlist(datalist)) == "type1or2")] == 2))
    numimm_type1 = length(datalist) - 1 - numimm_type2
    if(is.na(pars1[11]) == FALSE)
    {
      if(pars1[11] < numimm_type2/numimm | pars1[11] > (1 - numimm_type1 /numimm)) { return(-Inf) }
      datalist[[1]]$not_present_type2 = max(0,round(pars1[11] * numimm) - numimm_type2)
      datalist[[1]]$not_present_type1 = numimm - (length(datalist) - 1) - datalist[[1]]$not_present_type2
    }
    logp0_type1 = DAISIE_loglik_CS_choice(pars1 = pars1[1:5],pars2 = pars2,brts = datalist[[1]]$island_age,stac = 0,missnumspec = 0,methode = methode,CS_version = CS_version,abstolint = abstolint,reltolint = reltolint)
    logp0_type2 = DAISIE_loglik_CS_choice(pars1 = pars1[6:10],pars2 = pars2,brts = datalist[[1]]$island_age,stac = 0,missnumspec = 0,methode = methode,CS_version = CS_version,abstolint = abstolint,reltolint = reltolint)
    loglik = datalist[[1]]$not_present_type1 * logp0_type1 + datalist[[1]]$not_present_type2 * logp0_type2
    logcond = (cond == 1) * log(1 - exp((datalist[[1]]$not_present_type1 + numimm_type1) * logp0_type1 + (datalist[[1]]$not_present_type2 + numimm_type2) * logp0_type2))
  }
  loglik = loglik - logcond
  
  if(length(datalist) > 1)
  {
    for(i in 2:length(datalist))
    {
      if(datalist[[i]]$type1or2 == 1)
      {
        pars = pars1[1:endpars1]
      } else {
        pars = pars1[6:10]
      }
      loglik = loglik + DAISIE_loglik_CS_choice(pars1 = pars,pars2 = pars2,brts = datalist[[i]]$branching_times,stac = datalist[[i]]$stac,missnumspec = datalist[[i]]$missing_species,methode = methode,CS_version = CS_version,abstolint = abstolint,reltolint = reltolint)
    }
  }
  return(loglik)
}

DAISIE_integrate <- function(initprobs,tvec,rhs_func,pars,rtol,atol,method)
{
  if(length(pars) <= 7)
  {
    return(DAISIE_integrate_const(initprobs,tvec,rhs_func,pars,rtol,atol,method))
  } else {
    return(DAISIE_integrate_time(initprobs,tvec,rhs_func,pars,rtol,atol,method))
  }
}

DAISIE_integrate_const <- function(initprobs,tvec,rhs_func,pars,rtol,atol,method)
{
  if(as.character(body(rhs_func)[3]) == "lx <- (length(x) - 1)/2")
  {
    lx <- (length(initprobs) - 1)/2
    parsvec <- c(DAISIE_loglik_rhs_precomp(pars,lx))
    y <- DAISIE_ode_FORTRAN(initprobs,tvec,parsvec,atol,rtol,method,runmod = "daisie_runmod")
  } else if(as.character(body(rhs_func)[3]) == "lx <- (length(x))/3")
  {
    lx <- (length(initprobs))/3
    parsvec <- c(DAISIE_loglik_rhs_precomp(pars,lx))
    y <- DAISIE_ode_FORTRAN(initprobs,tvec,parsvec,atol,rtol,method,runmod = "daisie_runmod2")
  } else
  {
    stop('The integrand function is written incorrectly.')
  }
  return(y)
}

#' @useDynLib DAISIE
DAISIE_ode_FORTRAN <- function(
  initprobs,
  tvec,
  parsvec,
  atol,
  rtol,
  methode,
  runmod = "daisie_runmod"
)
{
  N <- length(initprobs)
  kk <- parsvec[length(parsvec)]
  if(runmod == "daisie_runmod")
  {
    lx <- (N - 1)/2
  } else if(runmod == "daisie_runmod2")
  {
    lx <- N/3
  }
  probs <- deSolve::ode(y = initprobs, parms = c(lx + 0.,kk + 0.), rpar = parsvec[-length(parsvec)], 
                        times = tvec, func = runmod, initfunc = "daisie_initmod", 
                        ynames = c("SV"), dimens = N + 2, nout = 1, outnames = c("Sum"), 
                        dllname = "DAISIE",atol = atol, rtol = rtol, method = methode)[,1:(N + 1)]
  return(probs)
}
