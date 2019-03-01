dec2bin = function(y,ly)
{
  stopifnot(length(y) == 1, mode(y) == 'numeric')
  q1 = (y / 2) %/% 1
  r = y - q1 * 2
  res = c(r)
  while(q1 >= 1)
  {
    q2 = (q1 / 2) %/% 1
    r = q1 - q2 * 2
    q1 = q2
    res = c(r, res)
  }
  res = c(rep(0,ly - length(res)),res)
  return(res)
}

dec2binmat = function(y)
{
  numrows = 2^y
  res = matrix(0,numrows,y)
  for(i in 0:(numrows-1))
  {
    res[i + 1,] = dec2bin(i,y)
  }
  return(res)
}

bin2dec <- function(y)
{
  res <- y %*% 2^((length(y) - 1):0)
  return(as.numeric(res))
}

kimat <- function(dec2binmatk)
{
  ki <- matrix(0,dim(dec2binmatk)[1],dim(dec2binmatk)[1])
  for(i in 2:dim(dec2binmatk)[1])
  {
    locationones <- which(dec2binmatk[i,] == 1)
    for(j in 1:length(locationones))
    {
      dec2binmatki <- dec2binmatk[i,]
      dec2binmatki[locationones[j]] = 0
      j2 <- 1 + bin2dec(dec2binmatki)
      ki[i,j2] <- 1 
    }
  }
  return(ki)
}

kmini = function(dec2binmatk,lxm1,lxm2,lxe,sysdim = dim(dec2binmatk)[1])
{
  posc = Matrix::rowSums(dec2binmatk)
  negc = log2(sysdim) - posc
  kmin = rep(negc,each = lxm1 * lxm2 * lxe)
  dim(kmin) = c(lxm1,lxm2,lxe,sysdim)
  ki = kimat(dec2binmatk)
  res = list(kmin = kmin,ki = ki)
  return(res)
}

nndivdep = function(lxm1,lxm2,lxe,sysdim,Kprime,M,k,l)
{
  nnl = c(0,0:(lxm1 + 1))
  nnm = c(0,0:(lxm2 + 1))
  nne = c(0,0,0:(lxe + 1))
  lnnl = length(nnl)
  lnnm = length(nnm)
  lnne = length(nne)
  nil2lxm1= 2:(lxm1 + 1)
  nil2lxm2 = 2:(lxm2 + 1)
  nil2lxe = 3:(lxe + 2)
  nn = rowSums(expand.grid(n1 = nnl,n2 = nnm,n3 = nne))
  dim(nn) = c(lnnl,lnnm,lnne)
  nn = replicate(sysdim,nn)
  nilm1 = rep(1,lxm1)
  nilm2 = rep(1,lxm2)
  nile = rep(1,lxe)
  allc = 1:sysdim
  divdepfac = pmax(array(0,dim = c(lxm1+3,lxm2+3,lxe+4,sysdim)),1 - (nn + k)/Kprime)
  divdepfacmin1 = pmax(array(0,dim = c(lxm1+3,lxm2+3,lxe+4,sysdim)),1 - (nn + k - 1)/Kprime)
  divdepfac = divdepfac[nil2lxm1,nil2lxm2,nil2lxe,allc]
  divdepfacmin1 = divdepfacmin1[nil2lxm1,nil2lxm2,nil2lxe,allc]
  Mminm = M - nn[nil2lxm1,nil2lxm2,nile,allc]
  lminm1 = l - nn[nil2lxm1,nilm2,nile,allc]
  Mminlminm2 = M - l - nn[nilm1,nil2lxm2,nile,allc]
  res = list(nn = nn,divdepfac = divdepfac,divdepfacmin1 = divdepfacmin1,Mminm = Mminm,lminm1 = lminm1,Mminlminm2 = Mminlminm2)
  return(res)
}

selectrows = function(sysdim,order)
{
  mat = NULL
  for(i in seq(1,sysdim/order,by = 2))
  {
    group1 = (i - 1) * order + (1:order)
    group2 = i * order + (1:order)
    mat = rbind(mat,cbind(group1,group2))
  }
  return(mat)
}

DAISIE_loglik_rhs_IW = function(t,x,pars)
{
  lac = pars[[1]][1]
  mu = pars[[1]][2]
  Kprime = pars[[1]][3]
  gam = pars[[1]][4]
  laa = pars[[1]][5]
  M = pars[[1]][6]
  k = pars[[2]]
  ddep = pars[[3]]
  lxm1 = pars[[4]]$lxm1
  lxm2 = pars[[4]]$lxm2
  lxe = pars[[4]]$lxe
  sysdim = pars[[4]]$sysdim
  kmin = pars[[5]]$kmin
  kplus = k - kmin
  ki = pars[[5]]$ki
  nn = pars[[6]]$nn
  divdepfac = pars[[6]]$divdepfac
  divdepfacmin1 = pars[[6]]$divdepfacmin1
  Mminm = pars[[6]]$Mminm
  Mminm[Mminm < 0] = 0
  lminm1minkminplus1 = pars[[6]]$lminm1 - kmin + 1
  lminm1minkminplus1[lminm1minkminplus1 < 0] = 0
  Mminlminm2plus1 = pars[[6]]$Mminlminm2 + 1
  Mminlminm2plus1[Mminlminm2plus1 < 0] = 0
  
  dim(x) = c(lxm1,lxm2,lxe,sysdim)
  xx = array(0,dim = c(lxm1+2,lxm2+2,lxe+3,sysdim))
  xx[2:(lxm1+1),2:(lxm2+1),3:(lxe+2),1:sysdim] = x
  nil2lxm1 = 2:(lxm1 + 1)
  nil2lxm2 = 2:(lxm2 + 1)
  nil2lxe = 3:(lxe + 2)
  allc = 1:sysdim
  nilm1 = rep(1,lxm1)
  nile = rep(1,lxe)
  nilm2 = rep(1,lxm2)
  
  if(sysdim == 1 & lxm1 > 1)
  {
    dim(Mminm) = c(lxm1,lxm2,lxe)
  } else 
  if(sysdim == 1 & lxm1 == 1)
  {
    dim(Mminm) = c(lxm2,lxe)
  } else
  if(sysdim > 1 & lxm1 == 1)
  {
    dim(Mminm) = c(lxm2,lxe,allc)
  }

  dx =  gam * divdepfacmin1 * lminm1minkminplus1 * xx[nil2lxm1-1,nil2lxm2,nil2lxe,allc] + #immigration
    gam * divdepfacmin1 * Mminlminm2plus1 * xx[nil2lxm1,nil2lxm2-1,nil2lxe,allc] + #immigration
    mu * nn[nil2lxm1+1,nilm2,nile,allc] * xx[nil2lxm1+1,nil2lxm2,nil2lxe,allc] + #extinction non-endemics
    mu * nn[nilm1,nil2lxm2+1,nile,allc] * xx[nil2lxm1,nil2lxm2+1,nil2lxe,allc] + #extinction non-endemics
    mu * nn[nilm1,nilm2,nil2lxe+1,allc] * xx[nil2lxm1,nil2lxm2,nil2lxe+1,allc] + #extinction endemics
    lac * divdepfacmin1 * nn[nil2lxm1+1,nilm2,nile,allc] * xx[nil2lxm1+1,nil2lxm2,nil2lxe-2,allc] + #cladogenesis non-endemics
    lac * divdepfacmin1 * nn[nilm1,nil2lxm2+1,nile,allc] * xx[nil2lxm1,nil2lxm2+1,nil2lxe-2,allc] + #cladogenesis non-endemics
    lac * divdepfacmin1 * nn[nilm1,nilm2,nil2lxe-1,allc] * xx[nil2lxm1,nil2lxm2,nil2lxe-1,allc] + #cladogenesis endemics
    2 * kplus * lac * divdepfacmin1 * xx[nil2lxm1,nil2lxm2,nil2lxe-1,allc] + #cladogenesis species in tree
    laa * nn[nil2lxm1+1,nilm2,nile,allc] * xx[nil2lxm1+1,nil2lxm2,nil2lxe-1,allc] + #anagenesis non-endemics
    laa * nn[nilm1,nil2lxm2+1,nile,allc] * xx[nil2lxm1,nil2lxm2+1,nil2lxe-1,allc] + #anagenesis non-endemics
    -(laa * (nn[nil2lxm1,nil2lxm2,nile,allc] + kmin) + (gam * divdepfac * Mminm) +
        (lac * divdepfac + mu) * (nn[nil2lxm1,nil2lxm2,nil2lxe,allc] + k)) * xx[nil2lxm1,nil2lxm2,nil2lxe,allc]
  if(sysdim > 1)
  {
    dx = dx +
      laa * tensor::tensor(xx[nil2lxm1,nil2lxm2,nil2lxe,allc],ki,4,2) + # anagenesis in colonizing lineage
      2 * lac * divdepfacmin1 * tensor::tensor(xx[nil2lxm1,nil2lxm2,nil2lxe-1,allc],ki,4,2) # cladogenesis in colonizing lineage
  }
  dim(dx) = c(sysdim * lxm1 * lxm2 * lxe,1)
  return(list(dx))
}



#' Computes the loglikelihood of the DAISIE model with island-wide
#' diversity-dependence given data and a set of model parameters
#' 
#' Computes the loglikelihood of the DAISIE model given colonization and
#' branching times for lineages on an island, and a set of model parameters for
#' the DAISIE model with island-wide diversity-dependence
#' 
#' The output is a loglikelihood value
#' 
#' @param pars1 Contains the model parameters: \cr \cr \code{pars1[1]}
#' corresponds to lambda^c (cladogenesis rate) \cr \code{pars1[2]} corresponds
#' to mu (extinction rate) \cr \code{pars1[3]} corresponds to K (clade-level
#' carrying capacity) \cr \code{pars1[4]} corresponds to gamma (immigration
#' rate) \cr \code{pars1[5]} corresponds to lambda^a (anagenesis rate) \cr
#' \code{pars1[6]} is optional; it may contain M, the total number of species
#' on the mainland \cr \cr
#' @param pars2 Contains the model settings \cr \cr \code{pars2[1]} corresponds
#' to lx = length of ODE variable x \cr \code{pars2[2]} corresponds to ddmodel
#' = diversity-dependent model, model of diversity-dependence, which can be one
#' of\cr \cr ddmodel = 0 : no diversity dependence \cr ddmodel = 1 : linear
#' dependence in speciation rate \cr ddmodel = 11: linear dependence in
#' speciation rate and in immigration rate \cr ddmodel = 2 : exponential
#' dependence in speciation rate\cr ddmodel = 21: exponential dependence in
#' speciation rate and in immigration rate\cr Only ddmodel = 11 is currently
#' implemented \cr \cr \code{pars2[3]} corresponds to cond = setting of
#' conditioning\cr \cr cond = 0 : conditioning on island age \cr cond = 1 :
#' conditioning on island age and non-extinction of the island biota \cr \cr
#' \code{pars2[4]} Specifies whether intermediate output should be provided,
#' because computation may take long. Default is 0, no output. A value of 1
#' means the parameters and loglikelihood are printed. A value of 2 means also
#' intermediate progress during loglikelihood computation is shown.
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
#' Default is "ode45"
#' @param abstolint Absolute tolerance of the integration
#' @param reltolint Relative tolerance of the integration
#' @return The loglikelihood
#' @author Rampal S. Etienne & Bart Haegeman
#' @seealso \code{\link{DAISIE_ML_IW}}, \code{\link{DAISIE_loglik_CS}},
#' \code{\link{DAISIE_sim}}
#' @references Valente, L.M., A.B. Phillimore and R.S. Etienne (2015).
#' Equilibrium and non-equilibrium dynamics simultaneously operate in the
#' Galapagos islands. Ecology Letters 18: 844-852.
#' @keywords models
#' @export DAISIE_loglik_IW
DAISIE_loglik_IW = function(
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
  
  if(is.na(pars2[4]))
  {
    pars2[4] = 0
  }
  if(is.null(datalist[[1]]$brts_table))
  {
    datalist = Add_brt_table(datalist)
  }
  brts = c(-abs(datalist[[1]]$brts_table[,1]),0)
  clade = datalist[[1]]$brts_table[,2]
  event = datalist[[1]]$brts_table[,3]
  pars1 = as.numeric(pars1)
  if(length(pars1) == 5)
  {
     np = datalist[[1]]$not_present
     if(is.null(np))
     {
        np = datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2
     } 
     if(is.null(np))
     {
        cat('Number of species not present is misspecified.\n')
        loglik = NA
        return(loglik)
     }
     M = length(datalist) - 1 + np
  } else
  if(length(pars1) == 6)
  {
    M = pars1[6] 
  } else
  {
    cat('pars1 should contain 5 or 6 values.\n')
    loglik = NA
    return(loglik)
  }
  
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
  l = 0

  if(min(pars1) < 0)
  {
    cat('One or more parameters are negative.\n')
    loglik = -Inf
    return(loglik)
  }
  
  endemic = 0
  nonendemic1 = 0
  nonendemic2 = 0
  for(i in 2:length(datalist))
  {
    endemic = endemic + (datalist[[i]]$stac == 5) 
    nonendemic1 = nonendemic1 + (datalist[[i]]$stac == 3)
    nonendemic2 = nonendemic2 + (datalist[[i]]$stac == 1)
  }
  
  if((ddep == 1 | ddep == 11) & (ceiling(Kprime) < nonendemic1 + nonendemic2 + endemic + length(brts) - 2))
  {
    cat("The value of K\' is incompatible with the number of species in the clade. Likelihood for this parameter set will be set to -Inf.\n")
    loglik = -Inf
    return(loglik)
  }
  if(ddep == 1 | ddep == 11)
  {
    lx = min(1 + ceiling(Kprime),DDD::roundn(pars2[1]) )
  } else {
    lx = DDD::roundn(pars2[1])
  }
  lxm1 = max(clade) + 1
  lxm2 = min(lx,M + 1)
  lxe = lx
  
  sysdimchange = 1
  sysdim = 1
  totdim = lxm1 * lxm2 * lxe * sysdim
  probs = rep(0,totdim)
  probs[1] = 1
  loglik = 0
  expandvec = NULL
  for(k in 0:(length(brts) - 2))
  {
    if(pars2[4] == 2)
    { 
      cat(paste('k = ',k,', sysdim = ',sysdim,'\n',sep = ''))
      utils::flush.console()
    }
    dime = list(lxm1 = lxm1,lxm2 = lxm2,lxe = lxe,sysdim = sysdim)
    if(sysdimchange == 1)
    {
      if(sysdim > 1)
      {
        dec2binmatk = dec2binmat(log2(sysdim))
        kmi = kmini(dec2binmatk,lxm1,lxm2,lxe,sysdim)
      } else if(sysdim == 1)
      {
        kmi = list(kmin = 0,ki = NULL)
      }
      sysdimchange = 0
    }
    nndd = nndivdep(lxm1,lxm2,lxe,sysdim,Kprime,M,k,l)
    parslist = list(pars = pars1,k = k,ddep = ddep,dime = dime,kmi = kmi,nndd = nndd)
    y = deSolve::ode(y = probs,times = brts[(k + 1):(k + 2)],func = DAISIE_loglik_rhs_IW,parms = parslist,rtol = reltolint,atol = abstolint,method = methode)
    probs = y[2,2:(totdim + 1)]
    cp = checkprobs2(NA,loglik,probs); loglik = cp[[1]]; probs = cp[[2]]      
    dim(probs) = c(lxm1,lxm2,lxe,sysdim)
    
    if(k < (length(brts) - 2))
    {
      divdepfac = nndd$divdepfac
      if(event[k + 2] == 1)
      {
        Mminlminm2 = nndd$Mminlminm2
        Mminlminm2[Mminlminm2 < 0] = 0
        probs = gam * divdepfac * Mminlminm2 * probs[,,,1:sysdim]
        probs = c(probs,rep(0,totdim))
        l = l + 1
        sysdim = sysdim * 2
        expandvec = c(expandvec,clade[k + 2])
        sysdimchange = 1
      } else 
      {
        probs = lac * divdepfac * probs[,,,1:sysdim]
        if(event[k + 2] == 2)
        {
          tocollapse = which(expandvec == clade[k + 2])
          sr = selectrows(sysdim,2^(tocollapse - 1))
          probs = probs[,,,sr[,1]] + probs[,,,sr[,2]]
          sysdim = sysdim / 2
          dim(probs) = c(lxm1,lxm2,lxe,sysdim)
          expandvec = expandvec[-tocollapse]
          sysdimchange = 1
        }
      }
      cp = checkprobs2(NA,loglik,probs); loglik = cp[[1]]; probs = cp[[2]]      
      totdim = lxm1 * lxm2 * lxe * sysdim
      dim(probs) = c(totdim,1)
    }
  }
  dim(probs) = c(lxm1,lxm2,lxe,sysdim)
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
  if(length(status) > 0)
  {
    decstatus = bin2dec(rev(status))
  } else
  {
    decstatus = 0
  }
  #print(status)
  numcol = length(datalist) - 1
  loglik = loglik - (lgamma(M + 1) - lgamma(M - numcol + 1) - lgamma(nonendemic2 + 1) - lgamma(endemic + 1))  
  #print(loglik + log(probs[1,,1,]))
  loglik = loglik + log(probs[1 + nonendemic1,1 + nonendemic2,1 + endemic,1 + decstatus])
  
  if(cond > 0)
  {
    sysdim = 1
    totdim = lxm1 * lxm2 * lxe * sysdim
    dime = list(lxm1 = lxm1,lxm2 = lxm2,lxe = lxe,sysdim = sysdim)
    probs = rep(0,totdim)
    probs[1] = 1
    kmi = list(kmin = 0,ki = NULL)
    nndd = nndivdep(lxm2,lxe,sysdim,Kprime,M,k = 0)
    parslist = list(pars = pars1,k = k,ddep = ddep,dime = dime,kmi = kmi,nndd = nndd)
    y = deSolve::ode(y = probs,times = brts[(k + 1):(k + 2)],func = DAISIE_loglik_rhs_IW,parms = parslist,rtol = reltolint,atol = abstolint,method = methode)
    probs = y[2,2:(totdim + 1)]
    dim(probs) = c(lxm1,lxm2,lxe,sysdim)
    logcond = log(1 - probs[1,1,1,1])
    loglik = loglik - logcond
  }
  if(pars2[4] >= 1)
  {
    s1 = sprintf('Parameters: %f %f %f %f %f',pars1[1],pars1[2],pars1[3],pars1[4],pars1[5])
    s2 = sprintf(', Loglikelihood: %f',loglik)
    cat(s1,s2,"\n",sep = "")
    utils::flush.console()
  }

  return(as.numeric(loglik))
}
