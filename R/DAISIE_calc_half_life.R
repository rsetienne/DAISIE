DAISIE_calc_half_life <- function(island_replicates, pars) {
  #calculate initial number of species for each simulation
  N0 <- list()
  for (i in 1:length(island_replicates)) {
    N0[[i]] <- sum(island_replicates$stt_tables[[i]][1, 2:4])
  }
  #Half way between initial species diversity and K
  spec_half <- list()
  for (i in 1:length(N0)) {
    spec_half[[i]] <- N0[[i]] - ((N0[[i]] - pars[3]) / 2)
    spec_half[[i]] <- round(spec_half[[i]], digits = 0)
  }
  #get total number of species through time
  total_stt <- list()
  for (i in 1:length(island_replicates)) {
    total_stt[[i]] <- apply(island_replicates[[i]]$stt_table[, 2:4], 1, sum)
  }
  #the row on which species richness equals half way between N0 and K
  row_spec_half <- list()
  for (i in 1:length(island_replicates)) {
    row_spec_half[[i]] <- min(which(total_stt[[i]] == spec_half[[i]]))
  }
  #the time take to reach the half-life
  half_life <- list()
  for (i in 1:length(island_replicates)) {
    half_life[[i]] <- as.numeric(island_replicates[[i]]$stt_table[1, 1]) -
      (island_replicates[[i]]$stt_table[[row_spec_half[[i]], 1]])
  }
  #if half-life has not been reached calculate using
  #Diamond 1972 exponential model
  if (length(row_t_half == 0)) {
    last_row <- nrow(island_replicates[[i]])
    species_at_present <- sum(island_replicates[[i]]$stt_table[last_row, 2:4])
    half_life <- 1 / -log( (species_at_present / pars[3]) / (N0 - pars[3]))
  }
  #calculate half-life from DAISIE_sim using a spline to compensate for sampling
  total_spec <- matrix(ncol = 1, nrow = nrow(spec_through_time_tables[[1]]))
  total_spec[, 1] <- sim[[1]][[1]]$stt_all[, 1]
  total_spec <- cbind(total_spec, sum_spec_through_time)
  #plot a spline for species diversity through time for each rep
  time_seq <- seq(0, max(total_spec[, 1]), 0.001)
  splines <- list()
  half_life_predict <- list()
  for (i in 2:ncol(total_spec)) {
    splines[[i]] <- smooth.spline(x = total_spec[, 1],
                                  y = total_spec[, i],
                                  df = 10)
    half_life_predict[[i]] <- predict(splines[[i]], time_seq)
    half_life_predict[[i]] <- cbind(half_life_predict[[i]]$x,
                                    half_life_predict[[i]]$y)
  }
  half_time <- matrix(ncol = 1, nrow = length(sim))
  for (i in 1:ncol(sum_spec_through_time)) {
    half_time[[i]] <- which.min(abs(sum_spec_through_time[, i] -
                                      spec_half[i, 1]))
    half_time[[i]] <- total_spec[half_time[i], 1]
  }
  #the time take to reach the half-life
  half_life <- matrix(nrow = length(sim), ncol = 1)
  for (i in 1:length(sim)) {
    half_life[i, 1] <- max(total_spec[, 1]) - (half_time[i])
  }
  return(half_life)
}