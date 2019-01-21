DAISIE_probdist_rhs = function(t,x,m)
{
   x = pmax(x,0)
   #print(t)
   #flush.console()
   nx = sqrt(length(x))
   dim(x) = c(nx,nx)     
   xx = matrix(0,nx+3,nx+3)
   xx[3:(nx+2),3:(nx+2)] = x
   # 3 is where we start to count
   dx = m[[1]] * xx[3:(nx+2),2:(nx+1)] + m[[2]] * xx[3:(nx+2),4:(nx+3)] + m[[3]] * xx[4:(nx+3),3:(nx+2)] + m[[4]] * xx[2:(nx+1),4:(nx+3)] + m[[5]] * xx[1:(nx+0),4:(nx+3)] + m[[6]] * xx[2:(nx+1),3:(nx+2)] - m[[7]] * xx[3:(nx+2),3:(nx+2)]
   dim(dx) = c(nx^2,1)
   return(list(dx))
}

DAISIE_probdist = function(pars1,pars2,tvec,initEI = c(0,0),initprobs = NULL)
{
   lac = pars1[1]
   mu = pars1[2]
   ga = pars1[4]
   laa = pars1[5]
   lx = pars2[1]
   M = pars2[2]
   abstol = 1e-16
   reltol = 1e-10
   nx1 = rep(-2:lx,lx + 3)
   nx1 = nx1 * (nx1 >= 0)
   dim(nx1) = c(lx + 3,lx + 3)
   nx2 = t(nx1)
   m = list()
   m[[1]] = ga * (M - nx2[3:(lx + 2),2:(lx + 1)])  # I - 1
   m[[2]] = mu * nx2[3:(lx + 2),4:(lx + 3)]        # I + 1
   m[[3]] = mu * nx1[4:(lx + 3),3:(lx + 2)]        # E + 1
   m[[4]] = laa * nx2[3:(lx + 2),4:(lx + 3)]       # I + 1
   m[[5]] = lac * nx2[3:(lx + 2),4:(lx + 3)]       # I + 1
   m[[6]] = lac * nx1[2:(lx + 1),3:(lx + 2)]       # E - 1
   m[[7]] = (mu + lac) * nx1[3:(lx + 2),3:(lx + 2)] + (mu + laa + lac) * nx2[3:(lx + 2),3:(lx + 2)] + ga * (M - nx2[3:(lx + 2),3:(lx + 2)])                       # E, I, I
   if(!is.null(initprobs))
   {
      probs = initprobs
   } else
   {
      probs = matrix(0,lx,lx)
      probs[initEI[1] + 1,initEI[2] + 1] = 1 
   }
   dim(probs) = c(lx * lx,1)
   y = ode(probs,c(0,tvec),DAISIE_probdist_rhs,m,rtol = reltol,atol = abstol, method = "ode45")
   return(y)
}

DAISIE_convertprobdist = function(pb)
{
   out = list()
   dime = dim(pb)
   for(i in 1:dime[1])
   {
      pb2 = pb[i,2:dime[2]]
      d = sqrt(dime[2] - 1)
      dim(pb2) = c(d,d)
      out[[i]] = pb2
   }
   out[[i + 1]] = pb[,1]
   return(out)
}

DAISIE_margprobdist = function(pars1,pars2,tvec,initEI = c(0,0),initprobs = NULL,pb = NULL)
{
   if(is.null(pb))
   {
      pb = DAISIE_probdist(pars1,pars2,tvec,initEI,initprobs)
   }
   lx = pars2[1]
   pbE = matrix(nrow = length(tvec) + 1,ncol = lx)
   pbI = matrix(nrow = length(tvec) + 1,ncol = lx)
   pbN = matrix(nrow = length(tvec) + 1,ncol = 2 * lx - 1)
   for(i in 1:(length(tvec) + 1))
   {
      pbEI = pb[i,2:(lx * lx + 1)]
      dim(pbEI) = c(lx,lx)
      pbE[i,1:lx] = rowSums(pbEI)
      pbI[i,1:lx] = colSums(pbEI)
      pbN[i,1:(2 * lx - 1)] = antidiagSums(pbEI)
   }
   out = list(pbE,pbI,pbN)
   names(out) = c("pE","pI","pN")
   return(out)
}

DAISIE_numcol_dist = function(pars1,pars2,tvec)
{
   y = DAISIE_probdist(pars1,c(pars2[1],1),tvec)
   lx = pars2[1]
   probs00 = y[2:(length(tvec) + 1),2]
   probstp = y[2,2:(lx * lx + 1)]
   probseq = y[3,2:(lx * lx + 1)]
   dim(probstp) = c(lx,lx)
   dim(probseq) = c(lx,lx)
   ee = rep(0:(lx - 1),lx)
   dim(ee) = c(lx,lx)
   expEtpapprox = sum(ee * probstp)
   expEINtp = DAISIE_ExpEIN(tvec[1],pars1,1)
   cat('The total sum of the probabilities at the first time is',sum(probstp),'\n')
   cat('The approximation for the expected number of endemics is',expEtpapprox,'\n')
   cat('The true value for the expected number of endemics is',expEINtp[[1]],'\n')
   expEteqapprox = sum(ee * probseq)
   expEINteq = DAISIE_ExpEIN(Inf,pars1,1)
   cat('The total sum of the probabilities at the second time is',sum(probstp),'\n')
   cat('The approximation for the expected number of endemics is',expEteqapprox,'\n')
   cat('The true value for the expected number of endemics is',expEINteq[[1]],'\n')
   flush.console()
   M = pars2[2]
   if(!is.na(pars1[11]))
   {
       Mnonfinches = M - round(pars1[11] * M)
   } else {
       Mnonfinches = M   
   }       
   pC = dbinom(0:Mnonfinches,Mnonfinches,1 - probs00)
   expC = Mnonfinches * (1 - probs00)
   cat('The approximation for the expected number of colonizations is',expC,'\n')   
   out = list(pC,expC,expEINtp,expEtpapprox,expEINteq,expEteqapprox)
   names(out) = list("pC","expC","expEINtp","expEtpapprox","expEINteq","expEteqapprox")
   return(out)
}
                                  
DAISIE_numcol = function(pars1,pars2,tvec,initEI = NULL)
{
   lx = pars2[1]
   M = pars2[2]
   nC = length(initEI)
   lt = length(tvec)
   unique_initEI = unique(initEI)
   nuC = length(unique_initEI)
   if(!is.na(pars1[11]))
   {
      Mnonfinches = M - round(pars1[11] * M) - nC
   } else {
      Mnonfinches = M - nC  
   }
   pC = matrix(nrow = lt,ncol = Mnonfinches + nC + 1)
   y = DAISIE_probdist(pars1,c(pars2[1],1),tvec,initEI = c(0,0),initprobs = NULL)
   probs00 = y[2:(lt + 1),2]       
   expC = Mnonfinches * (1 - probs00)
   lpC = Mnonfinches + 1
   for(j in 1:lt)
   {
      pC[j,1:lpC] = dbinom(0:Mnonfinches,Mnonfinches,1 - probs00[j])
   }
   if(nuC > 0)
   {
      for(i in 1:nuC)
      {
         abund_initEI = 0
         for(j in 1:nC)
         {
            if(prod(initEI[[j]] == unique_initEI[[i]]))
            {
               abund_initEI = abund_initEI + 1
            }
         }
         y = DAISIE_probdist(pars1,c(pars2[1],1),tvec,initEI = unique_initEI[[i]],initprobs = NULL)
         probs00 = y[2:(lt + 1),2]
         expC = expC + abund_initEI * (1 - probs00)
         lpC = lpC + abund_initEI
         for(j in 1:lt)
         {
            pC[j,1:lpC] = DDD::conv(pC[j,1:(lpC - abund_initEI)],dbinom(0:abund_initEI,abund_initEI,1 - probs00[j]))
         }
      }
   } 
   names(expC) = tvec
   colnames(pC) = 0:(lpC - 1)
   rownames(pC) = tvec
   out = list(expC,pC)
   names(out) = c("expC","pC")
   return(out)
}

DAISIE_KLdist = function(pars1,pars2,tvec)
{
   y = DAISIE_probdist(pars1,pars2,tvec)
   lx = pars2[1]
   M = pars2[2]
   tp = tvec[1]
   teq = tvec[2]
   probstp = y[2,2:(lx * lx + 1)]
   probseq = y[3,2:(lx * lx + 1)]
   dim(probstp) = c(lx,lx)
   dim(probseq) = c(lx,lx)
   ee = rep(0:(lx - 1),lx)
   dim(ee) = c(lx,lx)
   expEapprox = sum(ee * probstp)
   print(sum(probstp))
   print(expEapprox)
   print(DAISIE_ExpEIN(teq,pars1,M)[[1]])
   expEapprox = sum(ee * probseq)
   print(sum(probseq))
   print(expEapprox)
   print(DAISIE_ExpEIN(Inf,pars1,M)[[1]])
   KLdist = sum(probseq * log(probseq/probstp))
   return(KLdist)
}