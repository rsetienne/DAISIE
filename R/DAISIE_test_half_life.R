#Test DAISIE_half_life spline
#Spline_vec is a vector e.g. c(2,5,10,100,500,1000,5000,10000)

DAISIE_test_half_life_spline <- function (time, M, pars, replicates,
                                          nonoceanic, divdepmodel = 'CS', 
                                          prop_type2_pool = NA, 
                                          replicates_apply_type2 = TRUE, 
                                          sample_freq = 25, plot_sims = TRUE,
                                          spline_df, plot_splines = TRUE)
{
  sim_reps <- DAISIE_sim(time,M,pars,replicates,nonoceanic,divdepmodel) 
    
  #initial number of species for each simulation
  N0 <- matrix(nrow = length(sim_reps), ncol = 1)
  for (i in 1:length(sim_reps))
  {
    N0[i,1] <- sum(sim_reps[[i]][[1]]$stt_all[1,2:4])
  }
  
  #Half way between initial species diversity and K
  spec_half <- matrix(nrow = length(sim_reps), ncol = 1)
  for (i in 1:length(N0))
  {
    spec_half[i,1] <- N0[i] - ((N0[i] - pars[3])/2)
  }
  
  #get total number of species through time
  spec_through_time_tables <- list()
  for (i in 1:length(sim_reps))
  {
    spec_through_time_tables[[i]] <- sim_reps[[i]][[1]]$stt_all[,2:4]
  }
  
  sum_spec_through_time <- matrix(ncol = length(sim_reps),nrow = nrow(spec_through_time_tables[[1]]))
  for (i in 1:length(sim_reps))
  {
    sum_spec_through_time[,i] <- apply(X = spec_through_time_tables[[i]], MARGIN = 1, FUN = sum)
  } 
  
  total_spec <- matrix(ncol = 1, nrow = nrow(spec_through_time_tables[[1]]))
  total_spec[,1] <- sim_reps[[1]][[1]]$stt_all[,1]
  total_spec <- cbind(total_spec,sum_spec_through_time) 
  
  #plot a spline for species diversity through time for each rep
  time_seq <- seq(0,max(total_spec[,1]),0.001)
  splines <- list()
  half_life_predict <- list()
  for (i in 2:ncol(total_spec))
  {
    for (j in spline_df)
    {
      for (rep in 1:length(spline_df))
      {
        splines[[rep]] <- smooth.spline(x = total_spec[,1],y = total_spec[,i],df=j)
        half_life_predict[[rep]] <- predict(splines[[rep]], time_seq)
        half_life_predict[[rep]] <- cbind(half_life_predict[[rep]]$x,half_life_predict[[rep]]$y)
      }
    }
  }
  if (plot_splines == TRUE)
    {
      for (rep in 1:length(spline_df))
      {
      plot(rev(half_life_predict[[rep]][,1]),half_life_predict[[rep]][,2],ylab = "Species Diversity", xlab = "Time")
      }
    }
}





#Function to test how changing the sampling frequency affects the accurancy of the half-life prediction
#sample_freq_vec is a vector e.g. c(10, 25, 50, 100, 1000, 5000, 10000)
DAISIE_test_half_life_sample_freq <- function (time,M,pars,replicates,nonoceanic,divdepmodel = 'CS',
                                               prop_type2_pool = NA,replicates_apply_type2 = TRUE,
                                               sample_freq,plot_sims = TRUE, sample_freq_vec)
{
  half_life_list <- list()
  for (rep in 1:length(sample_freq_vec))
  {
    sim_reps <- DAISIE_sim(time,M,pars,replicates,nonoceanic,divdepmodel,prop_type2_pool,
                           replicates_apply_type2,sample_freq=sample_freq_vec[[rep]],plot_sims)
    
    #initial number of species for each simulation
    N0 <- matrix(nrow = length(sim_reps), ncol = 1)
    for (i in 1:length(sim_reps))
    {
      N0[i,1] <- sum(sim_reps[[i]][[1]]$stt_all[1,2:4])
    }
    
    #Half way between initial species diversity and K
    spec_half <- matrix(nrow = length(sim_reps), ncol = 1)
    for (i in 1:length(N0))
    {
      spec_half[i,1] <- N0[i] - ((N0[i] - pars[3])/2)
    }
    
    #get total number of species through time
    spec_through_time_tables <- list()
    for (i in 1:length(sim_reps))
    {
      spec_through_time_tables[[i]] <- sim_reps[[i]][[1]]$stt_all[,2:4]
    }
    
    sum_spec_through_time <- matrix(ncol = length(sim_reps),nrow = nrow(spec_through_time_tables[[1]]))
    for (i in 1:length(sim_reps))
    {
      sum_spec_through_time[,i] <- apply(X = spec_through_time_tables[[i]], MARGIN = 1, FUN = sum)
    } 
    
    total_spec <- matrix(ncol = 1, nrow = nrow(spec_through_time_tables[[1]]))
    total_spec[,1] <- sim_reps[[1]][[1]]$stt_all[,1]
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
    
    half_time <- matrix(ncol=1,nrow=length(sim_reps))
    for (i in 1:ncol(sum_spec_through_time))
    {
      half_time[[i]] <- which.min(abs(sum_spec_through_time[,i] - spec_half[i,1]))
      half_time[[i]] <- total_spec[half_time[i],1]
    }  
    
    #the time take to reach the half-life
    half_life <- matrix(nrow = length(sim_reps), ncol = 1)
    for (i in 1:length(sim_reps))
    {
      half_life[i,1] <- max(total_spec[,1]) - (half_time[i])
    }
    print(half_life)
  }
}
