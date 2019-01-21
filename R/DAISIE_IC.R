DAISIE_IC <- function(datalist,initparsopt,idparsopt,parsfix,idparsfix,endmc = 1000,res = 100,cond = 0,ddmodel = 0)
{
  MLE_obs <- DAISIE_ML(
     datalist = datalist,
     initparsopt = initparsopt,
     idparsopt = idparsopt,
     parsfix = parsfix,
     idparsfix = idparsfix,
     idparsnoshift = 6:10,
     res = res,
     ddmodel = ddmodel,
     cond = cond,
     eqmodel = 0,
     x_E = 0.95,
     x_I = 0.98,
     tol = c(1e-04, 1e-05, 1e-07),
     maxiter = 1000 * round((1.25)^length(idparsopt)),
     methode = "lsodes",
     optimmethod = 'subplex'
     )
  sims <- DAISIE_sim(
     time = datalist$island_age,
     M = datalist$not_present,         #add the number of species that are present
     pars = MLE_obs[1:5],
     replicates = endmc,
     format = TRUE,
     sample_freq = 1,
     plot_sims = FALSE
     )
  MLE <- rep(0,endmc)
  LL <- rep(0,endmc)
  for(mc in 1:endmc)
  {
    MLE[[mc]] <- DAISIE_ML(
       datalist = sims[[mc]],
       initparsopt = MLE_obs[idparsopt],
       idparsopt = idparsopt,
       parsfix = MLE_obs[parsfix],
       idparsfix = idparsfix,
       idparsnoshift = 6:10,
       res = res,
       ddmodel = ddmodel,
       cond = cond,
       eqmodel = 0,
       x_E = 0.95,
       x_I = 0.98,
       tol = c(1e-04, 1e-05, 1e-07),
       maxiter = 1000 * round((1.25)^length(idparsopt)),
       methode = "lsodes",
       optimmethod = 'subplex'
       )$loglik
    LL[[mc]] <- DAISIE_loglik_all(
       pars1 = MLE[[mc]][1:5],
       pars2 = c(res,ddmodel,cond,0),
       datalist = datalist,
       methode = "lsodes"
       )
  }
  WIC <- -2 * MLE_obs$loglik - 2 * mean(LL - MLE)
  AICb <- -2 * MLE_obs$loglik - 4 * mean(LL - MLE_obs)
  out <- list(WIC = WIC,AICb = AICb)
  return(out)
}