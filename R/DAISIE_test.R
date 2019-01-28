#' @title Tests DAISIE for data sets included in the paper by Valente et al
#' @description Tests DAISIE for data sets included in the paper by Valente et al.
#' @keywords model island biogeography
#' @export

DAISIE_test <- function()
{
  Galapagos_datalist_2types = NULL; rm(Galapagos_datalist_2types)
  Macaronesia_datalist = NULL; rm(Macaronesia_datalist)
  utils::data(Galapagos_datalist_2types, package = "DAISIE")
  pars1 = c(0.195442017,0.087959583,Inf,0.002247364,0.873605049,3755.202241,8.909285094,14.99999923,0.002247364,0.873605049,0.163)
  pars2 = c(100,11,0,0)
  loglik = DAISIE_loglik_all(pars1,pars2,Galapagos_datalist_2types)
  testthat::expect_equal(loglik,-61.7094829913735978)

  utils::data(Macaronesia_datalist, package = "DAISIE")
  background = c(0,1.053151832,Inf,0.052148979,0.512939011)
  Canaries = c(0.133766934,1.053151832,Inf,0.152763179,0.512939011)
  pars1 = rbind(background,Canaries,background, background)
  pars2 = c(100,0,0,0)
  loglik = 0
  for(i in 1:length(Macaronesia_datalist))
  {
    loglik = loglik + DAISIE_loglik_all(pars1[i,],pars2,Macaronesia_datalist[[i]],methode = "lsodes")
  }
  testthat::expect_equal(loglik,-454.9347833283220552)
}
