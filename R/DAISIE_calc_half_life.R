calc_half_life <- function(stt_tables) {
  #initial number of species for each simulation
  N0 <- matrix(nrow = length(stt_tables), ncol = 1)
  for (i in 1:length(stt_tables)) {
    N0[i, 1] <- sum(stt_tables[[i]][1, 2:4])
  }
  #Half way between initial species diversity and K
  spec_half <- matrix(nrow = length(stt_tables), ncol = 1)
  for (i in 1:length(N0)) {
    spec_half[i, 1] <- N0[i] - ( (N0[i] - pars[3]) / 2)
    }
  #get total number of species through time
  total_spec_tables <- list()
  total_spec <- list()
  for (i in 1:length(stt_tables)) {
    total_spec_tables[[i]] <- stt_tables[[i]][, 2:4]
    total_spec[[i]] <- apply(X = total_spec_tables[[i]], MARGIN = 1, FUN = sum)
  }
  for (i in 1:length(stt_tables)) {
    row_t_half <- list()
    row_t_half[i] <- min(which(total_spec[[i]] == spec_half[i]))
  }
  #if half-life has not been reached calculate using 
  #Diamond 1972 exponential model
  if (length(row_t_half == 0)) {
    last_row <- nrow(stt_tables)
    species_at_present <- sum(stt_tables[last_row, 2:4])
    half_life <- 1 / -log((species_at_present / pars[3]) / (N0 - pars[3]))
  }
  #the time take to reach the half-life
  half_life <- matrix(nrow = length(stt_tables), ncol = 1)
  for (i in 1:length(stt_tables)) {
    half_life[i, 1] <- totaltime - (stt_tables[[row_t_half, 1]])
  }
  return(half_life)
}