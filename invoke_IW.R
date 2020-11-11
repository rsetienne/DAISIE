library(DAISIE)
library(Rcpp)
library(RcppEigen)

IW <- function(methode = 'odeint::runge_kutta_fehlberg78') {
  utils::data(frogs_datalist, package = "DAISIE")
  pars1 = c(0.2, 0.1, 1000.1, 0.001, 0.3)
  pars2 = c(40, 11, 0, 0)
  loglik_IW = DAISIE::DAISIE_loglik_IW(
    pars1 = pars1,
    pars2 = pars2,
    datalist = frogs_datalist,
    methode = methode,
    abstolint = 1E-12,
    reltolint = 1E-10,
  )
  print(loglik_IW, digits=17)
  return(loglik_IW)
}

