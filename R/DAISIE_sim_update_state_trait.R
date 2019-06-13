#' Updates state of island given sampled event
#' 
#' Makes the event happen by updating island species matrix and species IDs.
#' What event happens is determined by the sampling in the algorithm.
#' @param timeval current time of simulation
#' @param totaltime simulated amount of time
#' @param possible_event numeric indicating what event will happen.
#' @param maxspecID current species IDs
#' @param mainland_spec number of mainland species with trait state 0
#' @param mainland_spec0 number of mainland species with trait state 1
#' @param island_spec A matrix with species on island (state of system at each time point)
#' @param stt_table A species-through-time table
#' @seealso \link{DAISIE_sim_core}
DAISIE_sim_update_state <- function(timeval,
                                    totaltime,
                                    possible_event,
                                    maxspecID,
                                    mainland_spec0,
                                    mainland_spec1,
                                    island_spec,
                                    stt_table)
{  
  # if (possible_event > 4) {
  #   # Nothing happens
  # }
  
  ##########################################
  #IMMIGRATION with trait0
  if(possible_event == 1)
  {  	
    
    
    colonist = DDD::sample2(mainland_spec0,1)      
    
    
    
    if(length(island_spec[,1]) != 0)
    {
      isitthere = which(island_spec[,1] == colonist)   
    } else
    {
      isitthere = c()  
    }
    
    if(length(isitthere) == 0)   
    {
      
      island_spec = rbind(island_spec,c(colonist,0,colonist,timeval,"I",NA,NA,NA))
      
      
      
    }
    if(length(isitthere) != 0)   
    {
      island_spec[isitthere,] = c(colonist,0,colonist,timeval,"I",NA,NA,NA)
      
    }
  }
  
  #######################
  ##immigration with trait 1      
  if(possible_event == 2)
  {  	
    colonist = DDD::sample2(mainland_spec1,1)    
    
    
    if(length(island_spec[,1]) != 0)
    {
      isitthere = which(island_spec[,1] == colonist)
    } else
    {
      isitthere = c()
    }
    
    if(length(isitthere) == 0)
    {
      island_spec = rbind(island_spec,c(colonist,1,colonist,timeval,"I",NA,NA,NA))   
    }
    
    if(length(isitthere) != 0)
    {
      island_spec[isitthere,] = c(colonist,1,colonist,timeval,"I",NA,NA,NA)    
    }
  }
  
  ##########################################        
  #EXTINCTION with trait 0 
  if(possible_event == 3)
  { 
    
    state0_spec = which(island_spec[,2] == "0")
    extinct = DDD::sample2(state0_spec,1)    
    
    typeofspecies = island_spec[extinct,5]   
    
    if(typeofspecies == "I")     
    {
      island_spec = island_spec[-extinct,]
    } else
      #remove immigrant
      
      if(typeofspecies == "A")
      {
        island_spec = island_spec[-extinct,]
      } else
        #remove anagenetic
        
        if(typeofspecies == "C")               
        {
          #remove cladogenetic
          
          #first to find species with same ancestor and arrival time          
          
          sisters = intersect(which(island_spec[,3] == island_spec[extinct,3]),which(island_spec[,4] == island_spec[extinct,4]))
          survivors = sisters[which(sisters != extinct)]   
          
          if(length(sisters) == 2)   
          {
            #survivors status becomes anagenetic	  
            island_spec[survivors,5] = "A"  
            island_spec[survivors,c(6,7)] = c(NA,NA)    
            island_spec[survivors,8] = "Clado_extinct"
            island_spec = island_spec[-extinct,]
          } else if(length(sisters) >= 3)
          {		
            numberofsplits = nchar(island_spec[extinct,6])     
            
            mostrecentspl = substring(island_spec[extinct,6],numberofsplits)  
            
            if(mostrecentspl == "B")
            { 
              sistermostrecentspl = "A"
            } else if(mostrecentspl == "A")
            {
              sistermostrecentspl = "B"
            }
            
            motiftofind = paste(substring(island_spec[extinct,6],1,numberofsplits-1),sistermostrecentspl,sep = "")
            
            possiblesister = survivors[which(substring(island_spec[survivors,6],1,numberofsplits) == motiftofind)]
            
            #different rules depending on whether a B or A is removed. B going extinct is simpler because it only 
            #carries a record of the most recent speciation			
            if(mostrecentspl == "A")
            {								
              #change the splitting date of the sister species so that it inherits the early splitting that used to belong to A.
              tochange = possiblesister[which(island_spec[possiblesister,7] == max(as.numeric(island_spec[possiblesister,7])))]
              island_spec[tochange,7] = island_spec[extinct,7]	 
            }
            
            #remove the offending A/B from these species
            island_spec[possiblesister,6] = paste(substring(island_spec[possiblesister,6],1,numberofsplits - 1),
                                                  substring(island_spec[possiblesister,6],numberofsplits + 1,
                                                            nchar(island_spec[possiblesister,6])),sep = "")	
            island_spec = island_spec[-extinct,]
          }
        }
    island_spec = rbind(island_spec)	
  }
  
  #############################
  
  #EXTINCTION with trait 1 
  if(possible_event == 4)
  { 
    #this chooses the row of species data to remove
    state1_spec = which(island_spec[,2] == "1")
    extinct = DDD::sample2(state1_spec,1)    
    
    typeofspecies = island_spec[extinct,5]   
    
    if(typeofspecies == "I")     
    {
      island_spec = island_spec[-extinct,]
    } else
      #remove immigrant
      
      if(typeofspecies == "A")
      {
        island_spec = island_spec[-extinct,]
      } else
        #remove anagenetic
        
        if(typeofspecies == "C")               
        {
          #remove cladogenetic
          
          #first to find species with same ancestor and arrival time          
          
          sisters = intersect(which(island_spec[,3] == island_spec[extinct,3]),which(island_spec[,4] == island_spec[extinct,4]))
          survivors = sisters[which(sisters != extinct)]   
          
          if(length(sisters) == 2)   
          {
            #survivors status becomes anagenetic	  
            island_spec[survivors,5] = "A"  
            island_spec[survivors,c(6,7)] = c(NA,NA)    
            island_spec[survivors,8] = "Clado_extinct"
            island_spec = island_spec[-extinct,]
          } else if(length(sisters) >= 3)
          {		
            numberofsplits = nchar(island_spec[extinct,6])     
            
            mostrecentspl = substring(island_spec[extinct,6],numberofsplits)  
            
            if(mostrecentspl == "B")
            { 
              sistermostrecentspl = "A"
            } else if(mostrecentspl == "A")
            {
              sistermostrecentspl = "B"
            }
            
            motiftofind = paste(substring(island_spec[extinct,6],1,numberofsplits-1),sistermostrecentspl,sep = "")
            
            possiblesister = survivors[which(substring(island_spec[survivors,6],1,numberofsplits) == motiftofind)]
            
            #different rules depending on whether a B or A is removed. B going extinct is simpler because it only 
            #carries a record of the most recent speciation			
            if(mostrecentspl == "A")
            {								
              #change the splitting date of the sister species so that it inherits the early splitting that used to belong to A.
              tochange = possiblesister[which(island_spec[possiblesister,7] == max(as.numeric(island_spec[possiblesister,7])))]
              island_spec[tochange,7] = island_spec[extinct,7]	  
            }
            
            #remove the offending A/B from these species
            island_spec[possiblesister,6] = paste(substring(island_spec[possiblesister,6],1,numberofsplits - 1),
                                                  substring(island_spec[possiblesister,6],numberofsplits + 1,
                                                            nchar(island_spec[possiblesister,6])),sep = "")	
            island_spec = island_spec[-extinct,]
          }
        }
    island_spec = rbind(island_spec)	
  }
  
  ##########################################
  #ANAGENESIS trait state  0
  if(possible_event == 5)
  {    
    immi_specs = intersect(which(island_spec[,2] == "0"), which(island_spec[,5] == "I"))
    
    #we only allow immigrants to undergo anagenesis
    if(length(immi_specs) == 1)
    {
      anagenesis = immi_specs
    } else if(length(immi_specs) > 1)
    {
      anagenesis = DDD::sample2(immi_specs,1)
    }
    
    maxspecID = maxspecID + 1                      
    island_spec[anagenesis,1] = maxspecID          
    island_spec[anagenesis,2] = "0"
    island_spec[anagenesis,5] = "A"                
    island_spec[anagenesis,8] = "Immig_parent"     
  }
  
  ######################
  
  ###anagenesis with trait1
  if(possible_event == 6)
  {    
    immi_specs = intersect(which(island_spec[,2] == "1"), which(island_spec[,5] == "I"))
    
    #we only allow immigrants to undergo anagenesis
    if(length(immi_specs) == 1)
    {
      anagenesis = immi_specs
    } else if(length(immi_specs) > 1)
    {
      anagenesis = DDD::sample2(immi_specs,1)
    }
    
    maxspecID = maxspecID + 1                      
    island_spec[anagenesis,1] = maxspecID          
    island_spec[anagenesis,2] = "1"
    island_spec[anagenesis,5] = "A"                
    island_spec[anagenesis,8] = "Immig_parent"     
  }
  
  
  
  
  
  ##########################################
  #CLADOGENESIS - this splits species into two new species - both of which receive    (with trait0)
  if(possible_event == 7)
  { 
    ##select a species with state 0 on the island
    state0_spec = which(island_spec[,2] == "0")
    
    tosplit = DDD::sample2(state0_spec,1)
    
    #if the species undergoes a cladogenesis
    if(island_spec[tosplit,5] == "C")
    {
      #for daughter A
      
      island_spec[tosplit,5] = "C"
      island_spec[tosplit,1] = maxspecID + 1  
      island_spec[tosplit,2] = "0"
      oldstatus = island_spec[tosplit,6]
      island_spec[tosplit,6] = paste(oldstatus,"A",sep = "")   
      #island_spec[tosplit,7] = timeval      
      island_spec[tosplit,8] = NA
      
      
      #for daughter B         
      island_spec = rbind(island_spec,c(maxspecID + 2,0,island_spec[tosplit,3],island_spec[tosplit,4],
                                        "C",paste(oldstatus,"B",sep = ""),timeval,NA))
      
      maxspecID = maxspecID + 2
    } else {
      #if the species that speciates is not cladogenetic
      
      #for daughter A
      
      island_spec[tosplit,5] = "C"
      island_spec[tosplit,1] = maxspecID + 1
      island_spec[tosplit,2] = "0"
      island_spec[tosplit,6] = "A"
      island_spec[tosplit,7] = island_spec[tosplit,4]   
      island_spec[tosplit,8] = NA
      
      
      #for daughter B
      island_spec = rbind(island_spec,c(maxspecID + 2,0,island_spec[tosplit,3],island_spec[tosplit,4],"C","B",timeval,NA))
      
      maxspecID = maxspecID + 2
    }
  } 
  ################
  #CLADOGENESIS - this splits species into two new species - both of which receive    (with trait0)
  if(possible_event == 8)
  { 
    ##select a species with state 0 on the island
    state1_spec = which(island_spec[,2] == "1")
    
    tosplit = DDD::sample2(state1_spec,1)
    
    #if the species undergoes a cladogenesis
    if(island_spec[tosplit,5] == "C")
    {
      #for daughter A
      
      island_spec[tosplit,5] = "C"
      island_spec[tosplit,1] = maxspecID + 1  
      island_spec[tosplit,2] = "1"
      oldstatus = island_spec[tosplit,6]
      island_spec[tosplit,6] = paste(oldstatus,"A",sep = "")   
      #island_spec[tosplit,7] = timeval      
      island_spec[tosplit,8] = NA
      
      
      #for daughter B         
      island_spec = rbind(island_spec,c(maxspecID + 2,1,island_spec[tosplit,3],island_spec[tosplit,4],
                                        "C",paste(oldstatus,"B",sep = ""),timeval,NA))
      
      maxspecID = maxspecID + 2
    } else {
      #if the species that speciates is not cladogenetic
      
      #for daughter A
      
      island_spec[tosplit,5] = "C"
      island_spec[tosplit,1] = maxspecID + 1
      island_spec[tosplit,2] = "1"
      island_spec[tosplit,6] = "A"
      island_spec[tosplit,7] = island_spec[tosplit,4]   
      island_spec[tosplit,8] = NA
      
      
      #for daughter B
      island_spec = rbind(island_spec,c(maxspecID + 2,1,island_spec[tosplit,3],island_spec[tosplit,4],"C","B",timeval,NA))
      
      maxspecID = maxspecID + 2
    }
  }
  
  ##########################################
  ##Transition 0 to 1
  if(possible_event == 9)
  {
    ## select a species with state 0
    state0_spec = which(island_spec[,2] == "0")
    totrans = DDD::sample2(state0_spec,1)
    
    island_spec[totrans,2] = "1"
    
  }
  
  
  ##########################################
  ##Transition 1 to 0
  if(possible_event == 10)
  {
    ## select a species with state 0
    state1_spec = which(island_spec[,2] == "1")
    totrans = DDD::sample2(state1_spec,1)
    
    island_spec[totrans,2] = "0"
    
  }  
  
  
 
  
  updated_state <- list(island_spec = island_spec, 
                        maxspecID = maxspecID, 
                        stt_table = stt_table)
  updated_state
}