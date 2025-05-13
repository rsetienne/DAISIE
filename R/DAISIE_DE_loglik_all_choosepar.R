DAISIE_DE_loglik_all_choosepar <- function(trparsopt,
                                           trparsfix,
                                           idparsopt,
                                           idparsfix,
                                           idparsnoshift,
                                           idparseq,
                                           pars2,
                                           datalist,
                                           methode,
                                           CS_version = list(model = 1, function_to_optimize = 'DAISIE'),
                                           abstolint = 1E-15,
                                           reltolint = 1E-15,
                                           equal_extinction = TRUE) {
  if(sum(idparsnoshift == (6:10)) != 5)
  {
    trpars1 <- rep(0,10)
  } else {
    trpars1 <- rep(0,5)
  }
  trpars1[idparsopt] <- trparsopt
  if(length(idparsfix) != 0)
  {
    trpars1[idparsfix] <- trparsfix
  }
  if(sum(idparsnoshift == (6:10)) != 5)
  {
    trpars1[idparsnoshift] <- trpars1[idparsnoshift - 5]
  }
  if(max(trpars1) > 1 | min(trpars1) < 0)
  {
    loglik <- -Inf
  } else {
    pars1 <- trpars1/(1 - trpars1)
    if (equal_extinction) {
      pars1[3] <- pars1[2]
      if(length(pars1) > 5) pars1[8] <- pars1[7]
    }
    if(min(pars1) < 0)
    {
      loglik <- -Inf
    } else {
      loglik <- DAISIE_DE_loglik_CS(pars1 = pars1,
                                    pars2 = pars2,
                                    datalist = datalist,
                                    methode = methode,
                                    abstolint,
                                    reltolint,
                                    equal_extinction)
    }
    if(is.nan(loglik) || is.na(loglik))
    {
      cat("There are parameter values used which cause numerical problems.\n")
      loglik <- -Inf
    }
  }
  return(loglik)
}



