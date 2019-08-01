#' Measures relaxation half-life from output 
#' from DAISIE_sim_core
#'
#' @param sim_core output from DAISIE_sim_core
#' @param pars parameter used to simulate data
#' 
#' @return a half life of the island
#' @export
#' @author Joshua Lambert
DAISIE_half_life <- function(sim_core, pars) {
  #initial number of species
  N0 <- sum(sim_core$stt_table[1, 2:4])
  #half-life of time taken to reach half way between initial species 
  #diversity and K
  N_half <- N0 - ((N0 - pars[3]) / 2)
  N_half <- round(N_half, digits = 0)
  #which row is the half-life number of species on
  stt_spec <- sim_core$stt_table[,2:4]
  sum_spec <- apply(X = test, MARGIN = 1, FUN = sum)
  row_t_half <- which(sum_spec == N_half)
  #if half-life has not been reached calculate
  if (length(row_t_half == 0)) {
    last_row <- nrow(sim_core$stt_table)
    species_at_present <- sum(sim_core$stt_table[last_row, 2:4])
    half_life = 1 / -log((species_at_present / pars[3]) / (N0 - pars[3]))
  }
  #what is the time when the half-life is reached
  time <- sim_core$stt_table[[1, 1]]
  #the time take to reach the half-life
  half_life <- time - (sim_core$stt_table[[t_half, 1]])
  return(half_life)
}

#' Measures relaxation average half-life from output 
#' from DAISIE_sim
#'
#' @param sim output from DAISIE_sim
#'
#' @return an average half life of the island
#' @export
#' @author Joshua Lambert
DAISIE_avg_half_life <- function (sim)
{
  #initial number of species for each simulation
  N0 <- matrix(nrow = length(sim), ncol = 1)
  for (i in 1:length(sim)) {
    N0[i, 1] <- sum(sim[[i]][[1]]$stt_all[1, 2:4])
  }
  
  #Half way between initial species diversity and K
  spec_half <- matrix(nrow = length(sim), ncol = 1)
  for (i in 1:length(N0)) {
    spec_half[i, 1] <- N0[i] - ((N0[i] - pars[3])/2)
  }
  
  #get total number of species through time
  spec_through_time_tables <- list()
  for (i in 1:length(sim))
  {
    spec_through_time_tables[[i]] <- sim[[i]][[1]]$stt_all[,2:4]
  }
  
  sum_spec_through_time <- matrix(ncol = length(sim),nrow = nrow(spec_through_time_tables[[1]]))
  for (i in 1:length(sim))
  {
    sum_spec_through_time[,i] <- apply(X = spec_through_time_tables[[i]], MARGIN = 1, FUN = sum)
  }
  
  total_spec <- matrix(ncol = 1, nrow = nrow(spec_through_time_tables[[1]]))
  total_spec[,1] <- sim[[1]][[1]]$stt_all[,1]
  total_spec <- cbind(total_spec,sum_spec_through_time)
  
  #plot a spline for species diversity through time for each rep
  time_seq <- seq(0,max(total_spec[,1]),0.001)
  splines <- list()
  half_life_predict <- list()
  for (i in 2:ncol(total_spec))
  {
    splines[[i]] <- smooth.spline(x = total_spec[,1],y = total_spec[,i],df=10)
    half_life_predict[[i]] <- predict(splines[[i]], time_seq)
    half_life_predict[[i]] <- cbind(half_life_predict[[i]]$x,half_life_predict[[i]]$y)
  }
  
  half_time <- matrix(ncol=1,nrow=length(sim))
  for (i in 1:ncol(sum_spec_through_time))
  {
    half_time[[i]] <- which.min(abs(sum_spec_through_time[,i] - spec_half[i,1]))
    half_time[[i]] <- total_spec[half_time[i],1]
  }
  
  #the time take to reach the half-life
  half_life <- matrix(nrow = length(sim), ncol = 1)
  for (i in 1:length(sim))
  {
    half_life[i,1] <- max(total_spec[,1]) - (half_time[i])
  }
  
  return(half_life)
}