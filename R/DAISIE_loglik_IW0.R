kmini0 = function(dec2binmatk,lxm,lxe,sysdim = dim(dec2binmatk)[1])
{
   posc = Matrix::rowSums(dec2binmatk)
   negc = log2(sysdim) - posc
   kmin = rep(negc,each = lxm * lxe)
   dim(kmin) = c(lxm,lxe,sysdim)
   ki = kimat(dec2binmatk)
   res = list(kmin = kmin,ki = ki)
   return(res)
}

nndivdep0 = function(lxm,lxe,sysdim,Kprime,M,k)
{
  nnm = c(0,0:(lxm + 1))
  nne = c(0,0,0:(lxe + 1))
  lnnm = length(nnm)
  lnne = length(nne)
  nn = matrix(nnm,nrow = lnnm,ncol = lnne,byrow = F) +
       matrix(nne,nrow = lnnm,ncol = lnne,byrow = T)
  nn = replicate(sysdim,nn)
  nil2lxm = 2:(lxm + 1)
  nil2lxe = 3:(lxe + 2)
  nile = rep(1,lxe)
  allc = 1:sysdim
  divdepfac = pmax(array(0,dim = c(lxm+3,lxe+4,sysdim)),1 - (nn + k)/Kprime)
  divdepfacmin1 = pmax(array(0,dim = c(lxm+3,lxe+4,sysdim)),1 - (nn + k - 1)/Kprime)
  divdepfac = divdepfac[nil2lxm,nil2lxe,allc]
  divdepfacmin1 = divdepfacmin1[nil2lxm,nil2lxe,allc]
  Mminm = M - nn[nil2lxm,nile,allc]
  res = list(nn = nn,divdepfac = divdepfac,divdepfacmin1 = divdepfacmin1,Mminm = Mminm)
  return(res)
}

DAISIE_loglik_rhs_IW0 = function(t,x,pars)
{
  lac = pars[[1]][1]
  mu = pars[[1]][2]
  Kprime = pars[[1]][3]
  gam = pars[[1]][4]
  laa = pars[[1]][5]
  M = pars[[1]][6]
  kk = pars[[2]]
  ddep = pars[[3]]
  lxm = pars[[4]]$lxm
  lxe = pars[[4]]$lxe
  sysdim = pars[[4]]$sysdim
  kmin = pars[[5]]$kmin
  kplus = kk - kmin
  ki = pars[[5]]$ki
  nn = pars[[6]]$nn
  divdepfac = pars[[6]]$divdepfac
  divdepfacmin1 = pars[[6]]$divdepfacmin1
  Mminm = pars[[6]]$Mminm

  dim(x) = c(lxm,lxe,sysdim)
  xx = array(0,dim = c(lxm+2,lxe+3,sysdim))
  xx[2:(lxm+1),3:(lxe+2),1:sysdim] = x
  nil2lxm = 2:(lxm + 1)
  nil2lxe = 3:(lxe + 2)
  allc = 1:sysdim
  nile = rep(1,lxe)
  nilm = rep(1,lxm)
  Mminm2 = Mminm + 1 - kmin
  Mminm2[Mminm2 < 0] = 0
  Mminm[Mminm < 0] = 0
  if(sysdim == 1)
  {
    dim(Mminm) = c(lxm,lxe)
  }
  dx =  gam * divdepfacmin1 * Mminm2 * xx[nil2lxm-1,nil2lxe,allc] + #immigration
        mu * nn[nil2lxm+1,nile,allc] * xx[nil2lxm+1,nil2lxe,allc] + #extinction non-endemics
        mu * nn[nilm,nil2lxe+1,allc] * xx[nil2lxm,nil2lxe+1,allc] + #extinction endemics
        lac * divdepfacmin1 * nn[nil2lxm+1,nile,allc] * xx[nil2lxm+1,nil2lxe-2,allc] + #cladogenesis non-endemics
        lac * divdepfacmin1 * nn[nilm,nil2lxe-1,allc] * xx[nil2lxm,nil2lxe-1,allc] + #cladogenesis endemics
        2 * kplus * lac * divdepfacmin1 * xx[nil2lxm,nil2lxe-1,allc] + #cladogenesis species in tree
        laa * nn[nil2lxm+1,nile,allc] * xx[nil2lxm+1,nil2lxe-1,allc] + #anagenesis non-endemics
       -(laa * (nn[nil2lxm,nile,allc] + kmin) + (gam * divdepfac * Mminm) +
        (lac * divdepfac + mu) * (nn[nil2lxm,nil2lxe,allc] + kk)) * xx[nil2lxm,nil2lxe,allc]
  if(sysdim > 1)
  {
      dx = dx +
      laa * tensor::tensor(xx[nil2lxm,nil2lxe,allc],ki,3,2) + # anagenesis in colonizing lineage
      2 * lac * divdepfacmin1 * tensor::tensor(xx[nil2lxm,nil2lxe-1,allc],ki,3,2) # cladogenesis in colonizing lineage
  }
  dim(dx) = c(sysdim * lxm * lxe,1)
  return(list(dx))
}


DAISIE_loglik_IW0 = function(
  pars1,
  pars2,
  datalist,
  methode = "ode45",
  abstolint = 1E-16,
  reltolint = 1E-14
  )
{
  # pars1 = model parameters
  # - pars1[1] = lac = (initial) cladogenesis rate
  # - pars1[2] = mu = extinction rate
  # - pars1[3] = K = maximum number of species possible in the clade
  # - pars1[4] = gam = (initial) immigration rate
  # - pars1[5] = laa = (initial) anagenesis rate
  # - pars1[6] = M = number of mainland species
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
  
  brts = c(-abs(datalist[[1]]$brts_table[,1]),0)
  clade = datalist[[1]]$brts_table[,2]
  event = datalist[[1]]$brts_table[,3]
  pars1 = as.numeric(pars1)

  ddep = pars2[2]
  cond = pars2[3]
  
  lac = pars1[1]
  mu = pars1[2]
  Kprime = pars1[3]
  if(ddep == 0)
  {
    Kprime = Inf
  }
  gam = pars1[4]
  laa = pars1[5]
  M = pars1[6]
  
  if(min(pars1) < 0)
  {
    cat('One or more parameters are negative.\n')
    loglik = -Inf
    return(loglik)
  }
  if((ddep == 1 | ddep == 11) & ceiling(Kprime) < length(brts))
  {
    cat('The value of K\' is incompatible with the number of species in the clade. Likelihood for this parameter set will be set to -Inf.\n')
    loglik = -Inf
    return(loglik)
  }
  if(ddep == 1 | ddep == 11)
  {
    lx = min(1 + ceiling(Kprime),DDD::roundn(pars2[1]) )
  } else {
    lx = DDD::roundn(pars2[1])
  }
  lxm = min(lx,M + 1)
  lxe = lx
  sysdimchange = 1
  sysdim = 1
  totdim = lxm * lxe * sysdim
  probs = rep(0,totdim)
  probs[1] = 1
  loglik = 0
  expandvec = NULL
  for(k in 0:(length(brts) - 2))
  {
    if(pars2[4] == 2)
    { 
       print(paste('k = ',k,', sysdim = ',sysdim,sep = ''))
       utils::flush.console()
    }
    dime = list(lxm = lxm,lxe = lxe,sysdim = sysdim)
    if(sysdimchange == 1)
    {
      if(sysdim > 1)
      {
        dec2binmatk = dec2binmat(log2(sysdim))
        kmi = kmini0(dec2binmatk,lxm,lxe,sysdim)
      } else if(sysdim == 1)
      {
        kmi = list(kmin = 0,ki = NULL)
      }
      sysdimchange = 0
    }
    nndd = nndivdep0(lxm,lxe,sysdim,Kprime,M,k)
    parslist = list(pars = pars1,kk = k,ddep = ddep,dime = dime,kmi = kmi,nndd = nndd)
    y = deSolve::ode(y = probs,times = brts[(k + 1):(k + 2)],func = DAISIE_loglik_rhs_IW,parms = parslist,rtol = reltolint,atol = abstolint,method = methode)
    probs = y[2,2:(totdim + 1)]
    cp = checkprobs2(NA,loglik,probs); loglik = cp[[1]]; probs = cp[[2]]      
    dim(probs) = c(lxm,lxe,sysdim)
    
    if(k < (length(brts) - 2))
    {
      divdepfac = nndd$divdepfac
      if(event[k + 2] == 1)
      {
        #Mminmminl = nndd$Mminm - kmi$kmin
        Mminmminl = nndd$Mminm - clade[k + 1]
        Mminmminl[Mminmminl < 0] = 0
        probs = gam * divdepfac * Mminmminl * probs[,,1:sysdim]
        probs = c(probs,rep(0,totdim))

        sysdim = sysdim * 2
        expandvec = c(expandvec,clade[k + 2])
        sysdimchange = 1
      } else 
      {
        probs = lac * divdepfac * probs[,,1:sysdim]
        if(event[k + 2] == 2)
        {
          tocollapse = which(expandvec == clade[k + 2])
          sr = selectrows(sysdim,2^(tocollapse - 1))
          probs = probs[,,sr[,1]] + probs[,,sr[,2]]
          sysdim = sysdim / 2
          dim(probs) = c(lxm,lxe,sysdim)
          expandvec = expandvec[-tocollapse]
          sysdimchange = 1
        }
      }
      cp = checkprobs2(NA,loglik,probs); loglik = cp[[1]]; probs = cp[[2]]      
      totdim = lxm * lxe * sysdim
      dim(probs) = c(totdim,1)
      #print(head(probs,n = 5))
    }
  }
  dim(probs) = c(lxm,lxe,sysdim)
  expandedclades = which(pracma::histc(clade,1:length(clade))$cnt == 1)
  status = rep(0,lexpandedclades <- length(expandedclades))
  if(lexpandedclades > 0)
  {
    for(i in lexpandedclades:1)
    {
      if(datalist[[1 + expandedclades[i]]]$stac == 2)
      {
        status[i] = 1
      }  
    }
  }
  endemic = 0
  nonendemic = 0
  for(i in 2:length(datalist))
  {
      endemic = endemic + (datalist[[i]]$stac == 5) 
      nonendemic = nonendemic + (datalist[[i]]$stac == 1) + (datalist[[i]]$stac == 3)
  }
  if(length(status) > 0)
  {
    decstatus = bin2dec(status)
  } else
  {
    decstatus = 0
  }
  print(loglik + log(probs))
  loglik = loglik + log(probs[1 + nonendemic,1 + endemic,1 + decstatus])

  if(cond > 0)
  {
    sysdim = 1
    totdim = lxm * lxe * sysdim
    dime = list(lxm = lxm,lxe = lxe,sysdim = sysdim)
    probs = rep(0,totdim)
    probs[1] = 1
    kmi = list(kmin = 0,ki = NULL)
    nndd = nndivdep0(lxm,lxe,sysdim,Kprime,M,k = 0)
    parslist = list(pars = pars1,kk = k,ddep = ddep,dime = dime,kmi = kmi,nndd = nndd)
    y = deSolve::ode(y = probs,times = brts[(k + 1):(k + 2)],func = DAISIE_loglik_rhs_IW0,parms = parslist,rtol = reltolint,atol = abstolint,method = methode)
    probs = y[2,2:(totdim + 1)]
    dim(probs) = c(lxm,lxe,sysdim)
    logcond = log(1 - probs[1,1,1])
    loglik = loglik - logcond
  }
  
  if(pars2[4] > 0)
  {
    s1 = sprintf('Parameters: %f %f %f %f %f %d',pars1[1],pars1[2],pars1[3],pars1[4],pars1[5],pars1[6])
    s2 = sprintf(', Loglikelihood: %f',loglik)
    cat(s1,s2,"\n",sep = "")
    utils::flush.console()
  }
  return(as.numeric(loglik))
}

DAISIE_loglik_IW_M1 <- function(
  pars1,
  pars2,
  brts,
  stac,
  missnumspec,
  methode = "ode45",
  abstolint = 1E-16,
  reltolint = 1E-14
  )
{
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
  
  if(is.na(pars2[4]))
  {
    pars2[4] = 0
  }
  ddep = pars2[2]
  cond = pars2[3]
  brts = c(-abs(brts),0)
  pars1 = as.numeric(pars1)
  lac = pars1[1]
  mu = pars1[2]
  Kprime = pars1[3]
  if(ddep == 0)
  {
    Kprime = Inf
  }
  gam = pars1[4]
  laa = pars1[5]
  pars1[6] = 1
  M = pars1[6]
  
  if(min(pars1) < 0)
  {
    cat('One or more parameters are negative.\n')
    loglik = -Inf
    return(loglik)
  }
  if(ddep == 1 | ddep == 11)
  {
    lx = min(1 + ceiling(Kprime),round(pars2[1]) )
  } else
  {
    lx = DDD::roundn(pars2[1])
  }
  lxm = min(lx,M + 1)
  lxe = lx
  sysdimchange = 1
  sysdim = 1
  totdim = lxm * lxe * sysdim
  probs = rep(0,totdim)
  probs[1] = 1
  loglik = 0
  expandvec = NULL
  for(k in 0:(length(brts) - 2))
  {
    if(pars2[4] == 2)
    { 
      print(paste('k = ',k,', sysdim = ',sysdim,sep = ''))
      utils::flush.console()
    }
    dime = list(lxm = lxm,lxe = lxe,sysdim = sysdim)
    if(sysdimchange == 1)
    {
      if(sysdim > 1)
      {
        dec2binmatk = dec2binmat(log2(sysdim))
        kmi = kmini0(dec2binmatk,lxm,lxe,sysdim)
      } else if(sysdim == 1)
      {
        kmi = list(kmin = 0,ki = NULL)
      }
      sysdimchange = 0
    }
    nndd = nndivdep0(lxm,lxe,sysdim,Kprime,M,k)
    parslist = list(pars = pars1,kk = k,ddep = ddep,dime = dime,kmi = kmi,nndd = nndd)
    y = deSolve::ode(y = probs,times = brts[(k + 1):(k + 2)],func = DAISIE_loglik_rhs_IW0,parms = parslist,rtol = reltolint,atol = abstolint,method = methode)
    probs = y[2,2:(totdim + 1)]
    cp = checkprobs2(NA,loglik,probs); loglik = cp[[1]]; probs = cp[[2]]      
    dim(probs) = c(lxm,lxe,sysdim)
    
    if(k < (length(brts) - 2))
    {
      divdepfac = nndd$divdepfac
      if(k == 0)
      {
        #Mminmminl = nndd$Mminm - kmi$kmin
        Mminmminl = nndd$Mminm
        Mminmminl[Mminmminl < 0] = 0
        probs = gam * divdepfac * Mminmminl * probs[,,1:sysdim]
        probs = c(probs,rep(0,totdim))
        sysdim = 2
        sysdimchange = 1
      } else
      {
        probs = lac * divdepfac * probs[,,1:sysdim]
        if(k == 1)
        {
          probs = probs[,,1] + probs[,,2]
          sysdim = 1
          dim(probs) = c(lxm,lxe,sysdim)
          sysdimchange = 1
        }
      }
      cp = checkprobs2(NA,loglik,probs); loglik = cp[[1]]; probs = cp[[2]]      
      totdim = lxm * lxe * sysdim
      dim(probs) = c(totdim,1)
    }
  }
  dim(probs) = c(lxm,lxe,sysdim)
  endemic = (stac == 5)
  nonendemic = (stac == 1) + (stac == 3)
  decstatus = (stac == 2) * (sysdim > 1) #when stac = 4, decstatus = 0
  #print(probs)
  loglik = loglik + log(probs[1 + nonendemic,1 + endemic,1 + decstatus])
  
  if(pars2[4] > 0)
  {
    s1 = sprintf('Parameters: %f %f %f %f %f %d',pars1[1],pars1[2],pars1[3],pars1[4],pars1[5],pars1[6])
    s2 = sprintf(', Loglikelihood: %f',loglik)
    cat(s1,s2,"\n",sep = "")
    utils::flush.console()
  }
  return(as.numeric(loglik))
}

DAISIE_loglik_IW0_choosepar = function(
  trparsopt,
  trparsfix,
  idparsopt,
  idparsfix,
  pars2,
  datalist,
  methode,
  abstolint,
  reltolint
)
{
  trpars1 = rep(0,6)
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
      loglik = DAISIE_loglik_IW0(pars1,pars2,datalist,methode,abstolint,reltolint)
    }
    if(is.nan(loglik) || is.na(loglik))
    {
      cat("There are parameter values used which cause numerical problems.\n")
      loglik = -Inf
    }
  }
  return(loglik)
}