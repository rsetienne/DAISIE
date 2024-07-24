DAISIE_loglik_high_lambda <- function(pars1, brts, stac) {
   lbrts <- length(brts)
   if (brts[lbrts] == 0) {
     brts <- brts[-lbrts]
     lbrts <- length(brts)
   }
   gam <- pars1[4]
   if (stac == 0) {
     out <- -gam * brts[1]
     return(out)
   }
   if (stac == 2 | stac == 6) {
     N <- lbrts - 1 - (stac == 6)
     mu <- pars1[2]
     brtsdiff <- brts - c(brts[2:lbrts], 0)
     out <- -gam * brtsdiff[1] +
       log(gam) * (stac == 2) +
       log(1 - exp(-gam * brtsdiff[2])) * (stac == 6) +
       #log(-expm1(-gam * brtsdiff[2])) * (stac == 6) +
       #log1p(-exp(-gam * brtsdiff[2])) * (stac == 6) +
       log(N) +
       (N - 1) * log(mu) +
       lgamma(N) +
       - (N - 1) * log(N - 1) +
       - mu / (N - 1) * sum((1:N) * (0:(N - 1)) * brtsdiff[(2 + (stac == 6)):(N + 1 + (stac == 6))])
     return(out)
   } else return(-Inf) # stac == 1, 3, 4, 5, 7, 8, 9
   #Stac 6 to be checked still
   return(out)
}
