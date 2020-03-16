DAISIE_sim_update_state_trait_dependent <- function(timeval,
                                          totaltime,
                                          possible_event,
                                          maxspecID,
                                          mainland_spec,
                                          island_spec,
                                          stt_table,
                                          trait_pars)
{
  if (possible_event > 10) {
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
      island_spec = rbind(island_spec,c(colonist,colonist,timeval,"I",NA,NA,NA,1))
    }
    
    if (length(isitthere) != 0)
    {
      island_spec[isitthere,] = c(colonist,colonist,timeval,"I",NA,NA,NA,1)
    }
  }
  
  ##########################################
  #EXTINCTION
  if (possible_event == 2)
  {
    island_spec_state1 = which(island_spec[,8] == "1")
    extinct = DDD::sample2(island_spec_state1,1)
    
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
          tochange = possiblesister[which(island_spec[possiblesister,6] == min(as.numeric(island_spec[possiblesister,6])))]
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
    immi_specs = intersect(which(island_spec[,4] == "I"), which(island_spec[,8] == "1"))
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
    if(!is.null(trait_pars)){
      island_spec[anagenesis,8] = "1"
    }
  }
  
  ##########################################
  #CLADOGENESIS - this splits species into two new species - both of which receive
  if(possible_event == 4)
  {
    island_spec_state1 = which(island_spec[,8] == "1")
    tosplit = DDD::sample2(island_spec_state1,1)
    
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
      island_spec[tosplit,8] = "1"
      
      #for daughter B
      island_spec = rbind(island_spec,c(maxspecID + 2,island_spec[tosplit,2],island_spec[tosplit,3],
                                        "C",paste(oldstatus,"B",sep = ""),timeval,NA,1))
      maxspecID = maxspecID + 2
    } else {
      #if the species that speciates is not cladogenetic
      
      #for daughter A
      
      island_spec[tosplit,4] = "C"
      island_spec[tosplit,1] = maxspecID + 1
      island_spec[tosplit,5] = "A"
      island_spec[tosplit,6] = island_spec[tosplit,3]
      island_spec[tosplit,7] = NA
      island_spec[tosplit,8] = "1"
      
      #for daughter B
      island_spec = rbind(island_spec,c(maxspecID + 2,island_spec[tosplit,2],island_spec[tosplit,3],"C","B",timeval,NA,1))
      maxspecID = maxspecID + 2
    }
  }
  
  ##########################
  ##transition from state1 to state2
  if(possible_event == 5){
    ##select a species with trait state1
    island_spec_state1 = which(island_spec[,8] == "1")
    totrans = DDD::sample2(island_spec_state1,1)
    island_spec[totrans,8] = "2"
  }
  
  ##########################
  ##immigration with state2
  if (possible_event == 6)
  {
    mainland1 = length(mainland_spec)
    mainland2 = trait_pars$M2
    mainland_total = mainland1 + mainland2
    colonist = DDD::sample2((mainland1 + 1):mainland_total,1)
    
    if (length(island_spec[,1]) != 0)
    {
      isitthere = which(island_spec[,1] == colonist)
    } else
    {
      isitthere = c()
    }
    
    if (length(isitthere) == 0)
    {
      island_spec = rbind(island_spec,c(colonist,colonist,timeval,"I",NA,NA,NA,2))
    }
    
    if (length(isitthere) != 0)
    {
      island_spec[isitthere,] = c(colonist,colonist,timeval,"I",NA,NA,NA,2)
    }
  }
  
  ##########################################
  #EXTINCTION
  if (possible_event == 7)
  {
    island_spec_state2 = which(island_spec[,8] == "2")
    extinct = DDD::sample2(island_spec_state2,1)
    #this chooses the row of species data with state2 to remove
    
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
      sisters = intersect(which(island_spec[,2] == island_spec[extinct,2]), which(island_spec[,3] == island_spec[extinct,3]))
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
          tochange = possiblesister[which(island_spec[possiblesister,6] == min(as.numeric(island_spec[possiblesister,6])))]
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
  if(possible_event == 8)
  {
    immi_specs = intersect(which(island_spec[,4] == "I"), which(island_spec[,8] == "2"))
    
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
    if(!is.null(trait_pars)){
      island_spec[anagenesis,8] = "2"
    }
  }
  
  ##########################################
  #CLADOGENESIS - this splits species into two new species - both of which receive
  if(possible_event == 9)
  {
    
    island_spec_state1 = which(island_spec[,8] == "2")
    tosplit = DDD::sample2(island_spec_state1,1)
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
      if(!is.null(trait_pars)){
        island_spec[tosplit,8] = "2"
      }
      #for daughter B
      if(is.null(trait_pars)){
        island_spec = rbind(island_spec,c(maxspecID + 2,island_spec[tosplit,2],island_spec[tosplit,3],
                                          "C",paste(oldstatus,"B",sep = ""),timeval,NA))
      }else{
        island_spec = rbind(island_spec,c(maxspecID + 2,island_spec[tosplit,2],island_spec[tosplit,3],
                                          "C",paste(oldstatus,"B",sep = ""),timeval,NA,2))
      }
      
      
      maxspecID = maxspecID + 2
    } else {
      #if the species that speciates is not cladogenetic
      
      #for daughter A
      
      island_spec[tosplit,4] = "C"
      island_spec[tosplit,1] = maxspecID + 1
      island_spec[tosplit,5] = "A"
      island_spec[tosplit,6] = island_spec[tosplit,3]
      island_spec[tosplit,7] = NA
      if(!is.null(trait_pars)){
        island_spec[tosplit,8] = "2"
      }
      #for daughter B
      if(is.null(trait_pars)){
        island_spec = rbind(island_spec,c(maxspecID + 2,island_spec[tosplit,2],island_spec[tosplit,3],"C","B",timeval,NA))
      }else{
        island_spec = rbind(island_spec,c(maxspecID + 2,island_spec[tosplit,2],island_spec[tosplit,3],"C","B",timeval,NA,2))
      }
      maxspecID = maxspecID + 2
    }
  }
  
  
  ##########################
  ##transition from state2 to state1
  if(possible_event == 10){
    ##select a species with trait state1
    island_spec_state1 = which(island_spec[,8] == "2")
    totrans = DDD::sample2(island_spec_state1,1)
    island_spec[totrans,8] = "1"
  }
  
  
  
  if (possible_event <= 10 && totaltime >= timeval) {
    stt_table <- rbind(stt_table,
                       c(totaltime - timeval,
                         length(intersect(which(island_spec[,4] == "I"),which(island_spec[,8] == "1"))),    #nI1
                         length(intersect(which(island_spec[,4] == "A"),which(island_spec[,8] == "1"))),    #nA1
                         length(intersect(which(island_spec[,4] == "C"),which(island_spec[,8] == "1"))),    #nC1
                         length(intersect(which(island_spec[,4] == "I"),which(island_spec[,8] == "2"))),    #nI2
                         length(intersect(which(island_spec[,4] == "A"),which(island_spec[,8] == "2"))),    #nA2
                         length(intersect(which(island_spec[,4] == "C"),which(island_spec[,8] == "2")))))   #nC2
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
