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
    
    ## this part can I use:  colonist = DDD::sample2(1:mainland_n0,1)  ???
    
    if(length(island_spec[,1]) != 0)
    {
      isitthere = which(island_spec[,1] == colonist)   ## to check if this species has already immigrate to the island,which means it is a re-immigration event
    } else
    {
      isitthere = c()  ##if the island is empty,just immigrate
    }
    
    if(length(isitthere) == 0)   ## means this is the first time to colonize on the island, just add a new line into the list
    {
      
      island_spec = rbind(island_spec,c(colonist,0,colonist,timeval,"I",NA,NA,NA))
      
      
      
      # if(colonist <= mainland_spec0 )    ##means the ancestor trait state of colonist is state 0
      #   {
      #   island_spec = rbind(island_spec,c(colonist,0,colonist,timeval,"I",NA,NA,NA)) ## it is the current trait state 0 in the twice column.
      # }
      # else if (colonist > mainland_spec0)     ##means the ancestor trait state of colonist is state 1
      # {
      #   island_spec = rbind(island_spec,c(colonist,1,colonist,timeval,"I",NA,NA,NA)) ## it is the current trait state 1 in the twice column.
      # }
    }
    if(length(isitthere) != 0)   ##means res-immigration happens, change the original state to the new one
    {
      island_spec[isitthere,] = c(colonist,0,colonist,timeval,"I",NA,NA,NA)
      
      # if(colonist <= mainland_spec0) 
      # {
      #   island_spec[isitthere,] = c(colonist,0,colonist,timeval,"I",NA,NA,NA)    
      # }
      #  
      # else if(colonist > mainland_spec0) 
      #  {
      #    island_spec[isitthere,] = c(colonist,1,colonist,timeval,"I",NA,NA,NA)   
      #  }
      
    }
  }
  
  
  ##########################################        
  #EXTINCTION with trait 0 
  if(possible_event == 2)
  { 
    #this chooses the row of species data to remove
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
        
        if(typeofspecies == "C")               #######？？？？？？？How to find sister clades when it has transform event
        {
          #remove cladogenetic
          
          #first to find species with same ancestor and arrival time          
          
          sisters = intersect(which(island_spec[,3] == island_spec[extinct,3]),which(island_spec[,4] == island_spec[extinct,4]))
          survivors = sisters[which(sisters != extinct)]   
          
          if(length(sisters) == 2)   ##means before extinction, it only happens cladogenesis for one time.
          {
            #survivors status becomes anagenetic	  ##all this parts change the number of columns
            island_spec[survivors,5] = "A"  
            island_spec[survivors,c(6,7)] = c(NA,NA)    
            island_spec[survivors,8] = "Clado_extinct"
            island_spec = island_spec[-extinct,]
          } else if(length(sisters) >= 3)
          {		
            numberofsplits = nchar(island_spec[extinct,6])     ### means how many ABs(times of branch)
            
            mostrecentspl = substring(island_spec[extinct,6],numberofsplits)  ##表示提取字符串的最后一个字符
            
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
              island_spec[tochange,7] = island_spec[extinct,7]	  #####??????????????????????????
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
  if(possible_event == 3)
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
    
    maxspecID = maxspecID + 1                      ##形成新物种后，更新物种库，所以总物种数增加1，并且新物种没有祖先状态，因此本身就是祖先I
    island_spec[anagenesis,1] = maxspecID          ##就是形成的新物种本身，但是祖先物种不发生改变
    island_spec[anagenesis,2] = "0"
    island_spec[anagenesis,5] = "A"                ##用A替换I，nI数量减少
    island_spec[anagenesis,8] = "Immig_parent"     ##来自于immigration
  }
  
  
  ##########################################
  #CLADOGENESIS - this splits species into two new species - both of which receive    (with trait0)
  if(possible_event == 4)
  { 
    ##select a species with state 0 on the island
    n0 = length(which(island_spec[,2] == "0"))
    
    tosplit = DDD::sample2(1:n0,1)
    
    #if the species undergoes a cladogenesis
    if(island_spec[tosplit,5] == "C")
    {
      #for daughter A
      
      island_spec[tosplit,5] = "C"
      island_spec[tosplit,1] = maxspecID + 1  
      island_spec[tosplit,2] = "0"
      oldstatus = island_spec[tosplit,6]
      island_spec[tosplit,6] = paste(oldstatus,"A",sep = "")   ###原来是A，现在变成AA和AB
      #island_spec[tosplit,7] = timeval      ###  branching time = current time
      island_spec[tosplit,8] = NA
      
      
      #for daughter B         here we assume that speciation will not cause state change.
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
      island_spec[tosplit,7] = island_spec[tosplit,4]   ###??????  branching time = clonisation time  ???
      island_spec[tosplit,8] = NA
      
      
      #for daughter B
      island_spec = rbind(island_spec,c(maxspecID + 2,0,island_spec[tosplit,3],island_spec[tosplit,4],"C","B",timeval,NA))
      
      maxspecID = maxspecID + 2
    }
  } 
  
  
  ##########################################
  ##Transition 0 to 1
  if(possible_event == 5)
  {
    ## select a species with state 0
    state0_spec = which(island_spec[,2] == "0")
    totrans = DDD::sample2(state0_spec,1)
    
    island_spec[totrans,2] = "1"
    
  }
  
  #######################
  ##immigration with trait 1      ##add a new event to general code
  if(possible_event == 6)
  {  	
    colonist = DDD::sample2(mainland_spec1,1)    ####change: mainland_sepc to mainland_spec1
    
    ##   colonist = DDD::sample2( (mainland_n0 +1) : mainland_ntotal,1) ???
    
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
  
  #############################
  
  #EXTINCTION with trait 1 
  if(possible_event == 7)
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
        
        if(typeofspecies == "C")               #######？？？？？？？How to find sister clades when it has transform event
        {
          #remove cladogenetic
          
          #first to find species with same ancestor and arrival time          
          
          sisters = intersect(which(island_spec[,3] == island_spec[extinct,3]),which(island_spec[,4] == island_spec[extinct,4]))
          survivors = sisters[which(sisters != extinct)]   
          
          if(length(sisters) == 2)   ##means before extinction, it only happens cladogenesis for one time.
          {
            #survivors status becomes anagenetic	  ##all this parts change the number of columns
            island_spec[survivors,5] = "A"  
            island_spec[survivors,c(6,7)] = c(NA,NA)    
            island_spec[survivors,8] = "Clado_extinct"
            island_spec = island_spec[-extinct,]
          } else if(length(sisters) >= 3)
          {		
            numberofsplits = nchar(island_spec[extinct,6])     ### means how many ABs(times of branch)
            
            mostrecentspl = substring(island_spec[extinct,6],numberofsplits)  ##表示提取字符串的最后一个字符
            
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
              island_spec[tochange,7] = island_spec[extinct,7]	  #####??????????????????????????
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
  
  ######################
  
  ###anagenesis with trait1
  if(possible_event == 8)
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
    
    maxspecID = maxspecID + 1                      ##形成新物种后，更新物种库，所以总物种数增加1，并且新物种没有祖先状态，因此本身就是祖先I
    island_spec[anagenesis,1] = maxspecID          ##就是形成的新物种本身，但是祖先物种不发生改变
    island_spec[anagenesis,2] = "1"
    island_spec[anagenesis,5] = "A"                ##用A替换I，nI数量减少
    island_spec[anagenesis,8] = "Immig_parent"     ##来自于immigration
  }
  
  ################
  #CLADOGENESIS - this splits species into two new species - both of which receive    (with trait0)
  if(possible_event == 9)
  { 
    ##select a species with state 0 on the island
    n1 = length(which(island_spec[,2] == "1")) 
    
    tosplit = DDD::sample2(1:n1,1)
    
    #if the species undergoes a cladogenesis
    if(island_spec[tosplit,5] == "C")
    {
      #for daughter A
      
      island_spec[tosplit,5] = "C"
      island_spec[tosplit,1] = maxspecID + 1  
      island_spec[tosplit,2] = "1"
      oldstatus = island_spec[tosplit,6]
      island_spec[tosplit,6] = paste(oldstatus,"A",sep = "")   ###原来是A，现在变成AA和AB
      #island_spec[tosplit,7] = timeval      ###  branching time = current time
      island_spec[tosplit,8] = NA
      
      
      #for daughter B         here we assume that speciation will not cause state change.
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