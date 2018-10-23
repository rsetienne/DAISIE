island_area <- function(t, totaltime, Apars, island_ontogeny)
{
  Tmax <- Apars$total_island_age
  Amax <- Apars$max_area
  Topt <- Apars$proportional_peak_t
  peak <- Apars$peak_sharpness
  proptime <- t/Tmax	
  # Constant
  if(is.null(island_ontogeny)) {
    return(Apars$max_area)
  }	
  if(island_ontogeny == "quadratic") {
    
    f <- Topt / (1 - Topt)
    a <- f * peak / (1 + f)
    b <- peak / (1 + f) 
    At <- Amax * proptime ^ a * (1 - proptime) ^ b / ((a / (a + b)) ^ a * (b / (a + b)) ^ b)
    return(At)}
}

DAISIE_loglik_rhs_time = function(t,x,parsvec)
{
  #parsvec[1:4] = Apars
  #parsvec[5] = lac0
  #parsvec[6:7] = mupars
  #parsvec[8] = K0
  #parsvec[9] = gam0
  #parsvec[10] = laa0
  #parsvec[11] = island_ontogeny
  #parsvec[12] = kk
  kk <- parsvec[length(parsvec)]
  lx <- (length(x) - 1)/2
  lnn <- lx + 4 + 2 * kk
  nn <- -2:(lx+2*kk+1)
  
  Apars <- parsvec[1:4]
  area <- island_area(t = t,Apars = Apars,island_ontogeny = parsvec[11])
  lacvec <- pmax(rep(0,lnn),parsvec[5] * (1 - nn/(area * parsvec[8])))
  X <- log(parsvec[6] / parsvec[7]) / log(0.1)
  mu <- parsvec[6] / ((area / Apars[1])^X)
  muvec <- mu * rep(1,lnn)
  gamvec <- pmax(rep(0,lnn),parsvec[9] * (1 - nn/(area * parsvec[8])))
  laavec <- parsvec[10] * rep(1,lnn)
  
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

DAISIE_loglik_rhs_time2 = function(t,x,parsvec)
{
  kk <- parsvec[length(parsvec)]
  lx <- (length(x) - 1)/3
  lnn <- lx + 4 + 2 * kk
  nn <- -2:(lx+2*kk+1)
  
  Apars <- parsvec[1:4]
  area <- island_area(t = t,Apars = Apars,island_ontogeny = parsvec[11])
  lacvec <- pmax(rep(0,lnn),parsvec[5] * (1 - nn/(area * parsvec[8])))
  X <- log(parsvec[6] / parsvec[7]) / log(0.1)
  mu <- parsvec[6] / ((area / Apars[1])^X)
  muvec <- mu * rep(1,lnn)
  gamvec <- pmax(rep(0,lnn),parsvec[9] * (1 - nn/(area * parsvec[8])))
  laavec <- parsvec[10] * rep(1,lnn)
  
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

divdepvec_time <- function(lacgam,pars1,lx,k1,ddep)
{
  if(lacgam == 'gam')
  {
    lacgam <- pars1[5]
  } else
  {
    lacgam <- pars1[9]
  }
  if(ddep == 1 | ddep == 11)
  {
    vec = pmax(rep(0,lx + 1),lacgam * (1 - ((0:lx)+k1) / K))
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

DAISIE_integrate_time <- function(initprobs,tvec,rhs_func,pars,rtol,atol,method)
{
  if(as.character(body(rhs_func)[3]) == "lx = (length(x) - 1)/2")
  {
    lx <- (length(initprobs) - 1)/2
    y <- ode(initprobs,tvec,parsvec,atol,rtol,method,runmod = "DAISIE_runmod")
  } else if(as.character(body(rhs_func)[3]) == "lx = (length(x))/3")
  {
    lx <- (length(initprobs))/3
    y <- ode(initprobs,tvec,parsvec,atol,rtol,method,runmod = "DAISIE_runmod2")
  } else
  {
    stop('The integrand function is written incorrectly.')
  }
  return(y)
}
