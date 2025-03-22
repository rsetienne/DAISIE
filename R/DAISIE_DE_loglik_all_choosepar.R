
###############################################################################
### function to calculate the log likelihood of all lineages in an island dataset
###############################################################################

DAISIE_DE_loglik_all_choosepar <- function(trparsopt,
                                           trparsfix,
                                           idparsopt,
                                           idparsfix,
                                           idparsnoshift,
                                           idparseq,
                                           pars2,
                                           datalist,
                                           methode,
                                           rtol = rtol,
                                           atol = atol,
                                           equal_extinction = TRUE) {
  if(sum(idparsnoshift == (6:10)) != 5)
  {
    trpars1 <- rep(0,10)
  } else {
    trpars1 <- rep(0,5)
  }
  trpars1[idparsopt] = trparsopt
  if(length(idparsfix) != 0)
  {
    trpars1[idparsfix] = trparsfix
  }
  if(sum(idparsnoshift == (6:10)) != 5)
  {
    trpars1[idparsnoshift] = trpars1[idparsnoshift - 5]
  }
  if(max(trpars1) > 1 | min(trpars1) < 0)
  {
    loglik <- -Inf
  } else {
    pars1 <- trpars1/(1 - trpars1)
    if (equal_extinction) {
      pars1[3] <- pars1[2]
    }
    if(min(pars1) < 0)
    {
      loglik <- -Inf
    } else {
      loglik <- DAISIE_DE_loglik_CS(pars1 = pars1,
                                    pars2 = pars2,
                                    datalist = datalist,
                                    methode = "lsodes",
                                    atol = atol,
                                    rtol = rtol)
    }
    if(is.nan(loglik) || is.na(loglik))
    {
      cat("There are parameter values used which cause numerical problems.\n")
      loglik <- -Inf
    }
  }
  return(loglik)
}
