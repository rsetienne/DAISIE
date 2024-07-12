DAISIE_loglik_high_lambda <- function(pars1, brts, stac) {
   lbrts <- length(brts)
   if (brts[lbrts] == 0) {
       brts <- brts[-lbrts]
       lbrts <- length(brts)
   }
   N <- lbrts - 1 - (stac == 6)
   mu <- pars1[2]
   gam <- pars1[4]
   brtsdiff <- brts - c(brts[2:(N + 1)], 0)
   if (stac == 0) {
      out <- -gam * brts[1]
   }
   if (stac == 2) {
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
   }
   if (stac == 7) {
     out <- log(1 - exp(-gamma * brtsdiff[1]))
   }
   if (stac == 1 | stac == 3 | stac == 4 | stac == 5 | stac >= 7) {
      out <- -Inf
   }
   #stac 6 still to be done.
   return(out)
}
