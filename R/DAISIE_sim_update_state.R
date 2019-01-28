#' Updates state of island given sampled event
#' 
#' Makes the event happen by updating island species matrix and species IDs.
#' What event happens is determined by the sampling in the algorithm.
#' @param timeval current time of simulation
#' @param totaltime simulated amount of time
#' @param possible_event numeric indicating what event will happen.
#' @param maxspecID current species IDs
#' @param mainland_spec number of mainland species
#' @param island_spec A matrix with species on island (state of system at each time point)
#' @param stt_table A species-through-time table
#' @seealso \link{DAISIE_sim_core}
DAISIE_sim_update_state <- function(timeval,
                                    totaltime,
                                    possible_event,
                                    maxspecID,
                                    mainland_spec,
                                    island_spec,
                                    stt_table)
{  
  if (possible_event > 4) {
    # Nothing happens
  }
  
  ##########################################
  #IMMIGRATION
  if (possible_event == 1)
  {
    colonist = DDD::sample2(mainland_spec,1)
    
    if (length(island_spec[,1]) != 0)
    {
      isitthere = which(island_spec[,1] == colonist)
    } else
    {
      isitthere = c()
    }
    
    if (length(isitthere) == 0)
    {
      island_spec = rbind(island_spec,c(colonist,colonist,timeval,"I",NA,NA,NA))
    }
    
    if (length(isitthere) != 0)
    {
      island_spec[isitthere,] = c(colonist,colonist,timeval,"I",NA,NA,NA)
    }
  }
  
  ##########################################
  #EXTINCTION
  if (possible_event == 2)
  { 	
    extinct = DDD::sample2(1:length(island_spec[,1]),1)
    #this chooses the row of species data to remove
    
    typeofspecies = island_spec[extinct,4]
    
    if(typeofspecies == "I")
    {
      island_spec = island_spec[-extinct,]
    }
    #remove immigrant
    
    if(typeofspecies == "A")
    {
      island_spec = island_spec[-extinct,]
    }
    #remove anagenetic
    
    if(typeofspecies == "C")
    {
      #remove cladogenetic
      #first find species with same ancestor AND arrival totaltime
      sisters = intersect(which(island_spec[,2] == island_spec[extinct,2]),which(island_spec[,3] == island_spec[extinct,3]))
      survivors = sisters[which(sisters != extinct)]
      
      if(length(sisters) == 2)
      {
        #survivors status becomes anagenetic	
        island_spec[survivors,4] = "A"
        island_spec[survivors,c(5,6)] = c(NA,NA)
        island_spec[survivors,7] = "Clado_extinct"
        island_spec = island_spec[-extinct,]
      }
      
      if(length(sisters) >= 3)
      {		
        numberofsplits = nchar(island_spec[extinct,5])
        
        mostrecentspl = substring(island_spec[extinct,5],numberofsplits)
        
        if(mostrecentspl=="B")
        { 
          sistermostrecentspl = "A"
        }
        if(mostrecentspl=="A")
        {
          sistermostrecentspl = "B"
        }
        
        motiftofind = paste(substring(island_spec[extinct,5],1,numberofsplits-1),sistermostrecentspl,sep = "")
        
        possiblesister = survivors[which(substring(island_spec[survivors,5],1,numberofsplits) == motiftofind)]
        
        #different rules depending on whether a B or A is removed. B going extinct is simpler because it only 
        #carries a record of the most recent speciation			
        if(mostrecentspl == "A")
        {								
          #change the splitting date of the sister species so that it inherits the early splitting that used to belong to A.
          tochange = possiblesister[which(island_spec[possiblesister,6] == max(as.numeric(island_spec[possiblesister,6])))]
          island_spec[tochange,6] = island_spec[extinct,6]	
        }
        
        #remove the offending A/B from these species
        island_spec[possiblesister,5] = paste(substring(island_spec[possiblesister,5],1,numberofsplits - 1),
                                              substring(island_spec[possiblesister,5],numberofsplits + 1,
                                                        nchar(island_spec[possiblesister,5])),sep = "")	
        island_spec = island_spec[-extinct,]
      }
    }
    island_spec = rbind(island_spec)	
  }
  
  ##########################################
  #ANAGENESIS
  if(possible_event == 3)
  {    
    immi_specs = which(island_spec[,4] == "I")
    
    #we only allow immigrants to undergo anagenesis
    if(length(immi_specs) == 1)
    {
      anagenesis = immi_specs
    }
    if(length(immi_specs) > 1)
    {
      anagenesis = DDD::sample2(immi_specs,1)
    }
    
    maxspecID = maxspecID + 1
    island_spec[anagenesis,4] = "A"
    island_spec[anagenesis,1] = maxspecID
    island_spec[anagenesis,7] = "Immig_parent"
  }
  
  ##########################################
  #CLADOGENESIS - this splits species into two new species - both of which receive 
  if(possible_event == 4)
  { 		
    tosplit = DDD::sample2(1:length(island_spec[,1]),1)
    
    #if the species that speciates is cladogenetic
    if(island_spec[tosplit,4] == "C")
    {
      #for daughter A
      
      island_spec[tosplit,4] = "C"
      island_spec[tosplit,1] = maxspecID + 1
      oldstatus = island_spec[tosplit,5]
      island_spec[tosplit,5] = paste(oldstatus,"A",sep = "")
      #island_spec[tosplit,6] = timeval
      island_spec[tosplit,7] = NA
      
      #for daughter B
      island_spec = rbind(island_spec,c(maxspecID + 2,island_spec[tosplit,2],island_spec[tosplit,3],
                                        "C",paste(oldstatus,"B",sep = ""),timeval,NA))
      
      maxspecID = maxspecID + 2
    } else {
      #if the species that speciates is not cladogenetic
      
      #for daughter A
      
      island_spec[tosplit,4] = "C"
      island_spec[tosplit,1] = maxspecID + 1
      island_spec[tosplit,5] = "A"
      island_spec[tosplit,6] = island_spec[tosplit,3]
      island_spec[tosplit,7] = NA
      
      #for daughter B
      island_spec = rbind(island_spec,c(maxspecID + 2,island_spec[tosplit,2],island_spec[tosplit,3],"C","B",timeval,NA))
      
      maxspecID = maxspecID + 2
    } 
  }
  
  
  if (possible_event <= 4 && totaltime >= timeval) {
    stt_table <- rbind(stt_table,
                       c(totaltime - timeval,
                         length(which(island_spec[,4] == "I")),
                         length(which(island_spec[,4] == "A")),
                         length(which(island_spec[,4] == "C"))))
  }
  
  updated_state <- list(island_spec = island_spec, 
                        maxspecID = maxspecID, 
                        stt_table = stt_table)
  updated_state
}

 
# !POTENTIALLY WRONG DUPLICATE FUNCTION! --------------------------------------- 
# DAISIE_ONEcolonist <- function(time,island_spec,stt_table)
# {
#   totaltime <- time
#   ### number of independent colonisations
#   uniquecolonisation <- as.numeric(unique(island_spec[,"Colonisation time (BP)"]))
#   number_colonisations <- length(uniquecolonisation) 
#   
#   ### if there is only one independent colonisation - anagenetic and cladogenetic
#   #species are classed as stac=2; immigrant classed as stac=4: 
#   if(number_colonisations == 1)
#   {
#     if(island_spec[1,"Species type"] == "I")
#     {
#       descendants <- list(stt_table = stt_table, 
#                           branching_times = c(totaltime,as.numeric(island_spec[1,"Colonisation time (BP)"])),
#                           stac = 4,
#                           missing_species = 0)
#     }
#     if(island_spec[1,"Species type"] == "A")
#     {
#       descendants <- list(stt_table = stt_table,
#                           branching_times = c(totaltime,as.numeric(island_spec[1,"Colonisation time (BP)"])),
#                           stac = 2,
#                           missing_species = 0)
#     } 
#     if(island_spec[1,"Species type"] == "C")
#     {
#       descendants <- list(stt_table = stt_table,
#                           branching_times = c(totaltime,rev(sort(as.numeric(island_spec[,"branching time (BP)"])))),
#                           stac = 2,
#                           missing_species = 0)
#     }
#   }
#   
#   ### if there are two or more independent colonisations, all species are classed as stac=3 and put within same list item: 
#   else if(number_colonisations > 1)
#   {
#     descendants <- list(stt_table = stt_table,
#                         branching_times = NA,stac = 2,missing_species = 0,
#                         other_clades_same_ancestor = list())
#     ### create table with information on other clades with same ancestor, but this information is not used in DAISIE_ML
#     oldest <- which(as.numeric(island_spec[,"Colonisation time (BP)"]) == max(as.numeric(island_spec[,"Colonisation time (BP)"])))
#     
#     oldest_table <- island_spec[oldest,]
#     if(class(oldest_table) == 'character')
#     { 
#       oldest_table <- t(as.matrix(oldest_table))
#     }
#     if(oldest_table[1,'Species type'] == 'A')
#     {
#       descendants$branching_times <- c(totaltime, as.numeric(oldest_table[1,"Colonisation time (BP)"]))
#     } else if(oldest_table[1,'Species type'] == 'C')
#     {
#       descendants$branching_times <- c(totaltime, rev(sort(as.numeric(oldest_table[,'branching time (BP)']))))
#     }
#     
#     youngest_table = island_spec[-oldest,]
#     if(class(youngest_table) == 'character')
#     {
#       youngest_table <- t(as.matrix(youngest_table))
#     }
#     
#     uniquecol <- as.numeric(unique(youngest_table[,"Colonisation time (BP)"]))
#     
#     descendants$missing_species <- length(which(youngest_table[,"Species type"]!='I'))
#     for(colonisation in 1:length(uniquecol))
#     {
#       descendants$other_clades_same_ancestor[[colonisation]] <- list(brts_miss = NA,species_type = NA)	
#       
#       samecolonisation <- which(as.numeric(youngest_table[,"Colonisation time (BP)"]) == uniquecol[colonisation])
#       
#       if(youngest_table[samecolonisation[1],"Species type"] == "I")
#       {
#         descendants$stac <- 3
#         descendants$other_clades_same_ancestor[[colonisation]]$brts_miss <- as.numeric(youngest_table[samecolonisation,"Colonisation time (BP)"])
#         descendants$other_clades_same_ancestor[[colonisation]]$species_type <- "I"
#       } else if(youngest_table[samecolonisation[1],"Species type"] == "A")
#       {
#         descendants$other_clades_same_ancestor[[colonisation]]$brts_miss <- as.numeric(youngest_table[samecolonisation,"Colonisation time (BP)"])
#         descendants$other_clades_same_ancestor[[colonisation]]$species_type <- "A"
#       } else if (youngest_table[samecolonisation[1],"Species type"] == "C")
#       {
#         descendants$other_clades_same_ancestor[[colonisation]]$brts_miss <- rev(sort(as.numeric(youngest_table[samecolonisation,"branching time (BP)"])))
#         descendants$other_clades_same_ancestor[[colonisation]]$species_type <- "C"
#       }
#     }
#   }
#   descendants
# }
