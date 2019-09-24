DAISIE_SR_sim_core <- function(time,mainland_n,pars,Tpars = NULL)
{
  totaltime <- time
  lac <- pars[1]
  mu <- pars[2]
  K <- pars[3]
  gam <- pars[4]
  laa <- pars[5]
  
  if(pars[4] == 0) 
  {
    stop('Rate of colonisation is zero. Island cannot be colonised.')
  }  
  
  timeval <- 0
  
  mainland_spec <- seq(1,mainland_n,1)
  maxspecID <- mainland_n
  
  island_spec <- c()
  stt_table <- matrix(ncol = 4)
  colnames(stt_table) <- c("Time","nI","nA","nC")
  stt_table[1,] <- c(totaltime,0,0,0)
  
  while(timeval < totaltime)
  {  	
    if(timeval < pars[11])
    { 
      lac <- pars[1]
      mu <- pars[2]
      K <- pars[3]
      gam <- pars[4]
      laa <- pars[5]
    } else
    { 
      lac <- pars[6]
      mu <- pars[7]
      K <- pars[8]
      gam <- pars[9]
      laa <- pars[10]
    }
    
    ext_rate <- mu * length(island_spec[,1])
    ana_rate <- laa * length(which(island_spec[,4] == "I"))
    clado_rate <- max(c(length(island_spec[,1]) * (lac * (1 -length(island_spec[,1])/K)),0),na.rm = T)
    immig_rate <- max(c(mainland_n * gam * (1 - length(island_spec[,1])/K),0),na.rm = T)
    
    totalrate <- ext_rate + clado_rate + ana_rate + immig_rate
    dt <- stats::rexp(1,totalrate)
    
    if ( timeval < pars[11] & ((timeval + dt) >= pars[11])  ) 
    {
      lac <- pars[6]
      mu <- pars[7]
      K <- pars[8]
      gam <- pars[9]
      laa <- pars[10]
      ext_rate <- mu * length(island_spec[,1])
      ana_rate <- laa * length(which(island_spec[,4] == "I"))
      clado_rate <- max(c(length(island_spec[,1]) * (lac * (1 -length(island_spec[,1])/K)),0),na.rm = T)
      immig_rate <- max( c(mainland_n * gam * (1 - length(island_spec[,1])/K), 0), na.rm = TRUE )
      totalrate <- ext_rate + clado_rate + ana_rate + immig_rate
      dt <- stats::rexp(1, totalrate)
      timeval <- pars[11] + dt
    } else 
    {
      timeval <- timeval + dt
    }
    
    possible_event <- sample(1:4,1,replace = FALSE,c(immig_rate,ext_rate,ana_rate,clado_rate))
    
    ##############
    if(timeval <= totaltime)
    { 
      new_state <- DAISIE_sim_update_state(timeval = timeval, 
                                           totaltime = totaltime,
                                           possible_event = possible_event,
                                           maxspecID = maxspecID,
                                           mainland_spec = mainland_spec,
                                           island_spec = island_spec,
                                           stt_table = stt_table,
                                           Tpars = Tpars)
      island_spec <- new_state$island_spec
      maxspecID <- new_state$maxspecID
    }
    stt_table <- rbind(stt_table,
                       c(totaltime - timeval,
                         length(which(island_spec[,4] == "I")),
                         length(which(island_spec[,4] == "A")),
                         length(which(island_spec[,4] == "C"))
                       )
    )
  }
  
  stt_table[nrow(stt_table),1] <- 0
  
  ############# 
  ### if there are no species on the island branching_times = island_age, stac = 0, missing_species = 0 
  if(length(island_spec[,1]) == 0)
  {
    island <- list(stt_table = stt_table, branching_times = totaltime, stac = 0, missing_species = 0)
  } else
  {
    cnames <- c("Species","Mainland Ancestor","Colonisation time (BP)",
                "Species type","branch_code","branching time (BP)","Anagenetic_origin")
    colnames(island_spec) <- cnames
    
    ### set ages as counting backwards from present
    island_spec[,"branching time (BP)"] <- totaltime - as.numeric(island_spec[,"branching time (BP)"])
    island_spec[,"Colonisation time (BP)"] <- totaltime - as.numeric(island_spec[,"Colonisation time (BP)"])
    
    if(mainland_n == 1)
    {
      island <- DAISIE_ONEcolonist(totaltime,island_spec,stt_table)
    } else if(mainland_n > 1)
    {  
      ### number of colonists present
      colonists_present <- sort(as.numeric(unique(island_spec[,'Mainland Ancestor'])))
      number_colonists_present <- length(colonists_present) 
      
      island_clades_info <- list()  
      for(i in 1:number_colonists_present)
      {
        subset_island <- island_spec[which(island_spec[,'Mainland Ancestor']==colonists_present[i]),] 
        if(class(subset_island) != 'matrix')
        {
          subset_island <- rbind(subset_island[1:7])
          colnames(subset_island) <- cnames
        }
        island_clades_info[[i]] <- DAISIE_ONEcolonist(totaltime,island_spec=subset_island,stt_table=NULL)
        island_clades_info[[i]]$stt_table <- NULL
      }
      island <- list(stt_table = stt_table, taxon_list = island_clades_info)
    }
  }
  return(island) 
}
