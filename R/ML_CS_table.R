ML_CS_store_params <-function(island_replicates,pars){
  lambda_c <- factor()
  mu <- factor()
  K <- factor()
  gamma <- factor()
  lambda_a <- factor()
  loglik <- factor()
  df <- factor()
  conv <- factor()
  DAISIEMLSUM_CS <- data.frame(lambda_c,mu,K,gamma,lambda_a
                               ,loglik,df,conv)
  for (i in 1:length(island_replicates)){
    DAISIEML_CS <- DAISIE::DAISIE_ML_CS(datalist=island_replicates[[i]], initparsopt = pars, ddmodel=0, idparsopt = 1:5, parsfix = NULL, idparsfix = NULL)
    DAISIEMLSUM_CS <- rbind(DAISIEMLSUM_CS, DAISIEML_CS)
  }
  return(DAISIEMLSUM_CS)
}

islands_sim_MLE <- function(DAISIEMLSUM_CS){
  islands_MLE = list()
  for(rep in 1:(length(DAISIEMLSUM_CS)-1)){
    lac <-DAISIEMLSUM_CS[rep,1]
    mu <- DAISIEMLSUM_CS[rep,2]
    k <- DAISIEMLSUM_CS[rep,3]
    gam <- DAISIEMLSUM_CS[rep,4]
    laa <- DAISIEMLSUM_CS[rep,5]
    pars = c(lac,mu,k,gam,laa)
    islands_MLE[[rep]] = DAISIE::DAISIE_sim(time=2,M=40,pars=pars,Tpars = NULL, replicates=10, divdepmodel = "CS")
  }
  samp_freq = length(islands_MLE[[rep]][[1]][[1]]$stt_all)
  sum_endemics = list()
  total = list()
  for(rep in 1:(length(DAISIEMLSUM_CS)-1)){
    sum_endemics[[rep]] <- islands_MLE[[rep]][[1]][[1]]$stt_all[samp_freq, "nA"] +
      islands_MLE[[rep]][[1]][[1]]$stt_all[samp_freq, "nC"]
    total[[rep]] <- islands_MLE[[rep]][[1]][[1]]$stt_all[samp_freq, "nA"] +
      islands_MLE[[rep]][[1]][[1]]$stt_all[samp_freq, "nC"] +
      islands_MLE[[rep]][[1]][[1]]$stt_all[samp_freq, "nI"]
  }
}
