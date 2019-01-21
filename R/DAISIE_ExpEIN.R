DAISIE_ExpEIN = function(t,pars,M,initEI = c(0,0))
{
   pars1 = pars
   lac = pars1[1]
   mu = pars1[2]
   ga = pars1[4]
   laa = pars1[5]
   if(!is.na(pars1[11]))
   {
       M2 = M - DDD::roundn(pars1[11] * M)
   } else {
       M2 = M
   }
   A = mu - lac
   B = lac + mu + ga + laa
   C = laa + 2 * lac + ga
   DD = laa + 2 * lac
   E0 = initEI[1]
   I0 = initEI[2] 
   if(t == Inf)
   {
      Imm = ga * M2 / B
      End = DD/A * Imm
   } else {
      #Imm = M2 * ga / B * (1 - exp(-B * t))
      #End = M2 * ga * (laa + 2 * lac) * (1/(A * B) - exp(-A*t) / (A * C) + exp(-B*t)/(B * C))
      Imm = M2 * ga / B - (M2 * ga / B - I0) * exp(-B * t)
      End = DD/C * (M2 * ga / A - M2 * ga/ B + (C / DD * E0 - M2 * ga / A + I0) * exp(-A * t) + (M2 * ga / B - I0) * exp(-B * t))
   }
   All = End + Imm
   expEIN = list(End,Imm,All)
   names(expEIN) = c("ExpE","ExpI","ExpN")
   return(expEIN)
}