#function to detect the time taken for the relaxation half life
#uses the object output from DAISIE_sim_core


#' Title
#'
#' @param time 
#' @param mainland_n 
#' @param pars 
#' @param nonoceanic 
#' @authorJosh Lambert 
#'
#' @return
#' @export
#'
#' @examples
DAISIE_sim_half_life <- function(time,mainland_n,pars,nonoceanic)
{
  #run simulation
  sim <- DAISIE_sim_core(time,mainland_n,pars,nonoceanic)
  
  #initial number of species
  N0 <- sum(sim$stt_table[1,2:4])
  
  #half-life of time taken to reach half way between initial species diversity and K
  t_half <- N0 - ((N0 - pars[3])/2)
  t_half <- round(t_half, digits = 0)
  
  #which row is the half-life number of species on
  test <- sim$stt_table[,2:4]
  row <- apply(X = test, MARGIN = 1, FUN = sum)
  row <- which(row == t_half)
  if (length(row) > 1)
  {
    row <- row[1]
  }
  
  #what is the time when the half-life is reached
  sim$stt_table[row,1]
  
  #the time take to reach the half-life
  half_life <- time - (sim$stt_table[row,1])
  
  return(half_life)
}

  
  
DAISIE_half_life <- function(sim)
{
  
  #initial number of species
  N0 <- sum(sim$stt_table[1,2:4])
  
  #half-life of time taken to reach half way between initial species diversity and K
  t_half <- N0 - ((N0 - pars[3])/2)
  t_half <- round(t_half, digits = 0)
  
  #which row is the half-life number of species on
  test <- sim$stt_table[,2:4]
  row <- apply(X = test, MARGIN = 1, FUN = sum)
  row <- which(row == t_half)
  if (length(row) > 1)
  {
    row <- row[1]
  }
  
  #what is the time when the half-life is reached
  sim$stt_table[row,1]
  
  #the time take to reach the half-life
  half_life <- time - (sim$stt_table[row,1])
  
  return(half_life)
}
 

DAISIE_avg_half_life <- function (sim_reps)
{
  #initial number of species for each simulation
  N0 <- matrix(nrow = length(sim_reps), ncol = 1)
  for (i in 1:length(sim_reps))
  {
    N0[i,1] <- sum(sim_reps[[i]][[1]]$stt_all[1,2:4])
  }
  
  #Half way between initial species diversity and K
  t_half <- matrix(nrow = length(sim_reps), ncol = 1)
  for (i in 1:length(N0))
  {
    t_half[i,1] <- N0[i] - ((N0[i] - pars[3])/2)
    t_half[i,1] <- round(t_half[i,1],digits = 0)
  }
  
  #which row is the half-life number of species on
  tables <- list()
  for (i in 1:length(sim_reps))
  {
    tables[[i]] <- sim_reps[[i]][[1]]$stt_all[,2:4]
  }
  
  rows <- list()
  for (i in 1:length(sim_reps))
  {
    rows[[i]] <- apply(X = tables[[i]], MARGIN = 1, FUN = sum)
  }  
  
  half_time_matrix <- matrix(nrow = length(sim_reps), ncol=1)
  for (i in length(sim_reps))
  {
    half_time <- which(rows[[i]] == t_half[i,1])
    if (length(half_time) > 1)
    {
      half_time <- half_time[1]
    }
    half_time_matrix[i,1]<-half_time
  }
  
  #the time take to reach the half-life
  half_life <- matrix(nrow = length(sim_reps), ncol = 1)
  for (i in 1:length(sim_reps))
  {
    half_life[i,1] <- time - (sim_reps[[i]][[1]]$stt_all[half_time_matrix[i],1])
  }
  
  return(half_life)
}



 
