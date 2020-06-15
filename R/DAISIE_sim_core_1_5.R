# This function produces an extra element per replicate when the
# island is empty at time = 0. Functionally this has no effect on the
# simulations, but care should be taken if using the length of objects to count
# the number of species present on the island.
DAISIE_sim_core_1_5 <- function(time,mainland_n,pars)
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

  island_spec = c()
  stt_table <- matrix(ncol = 4)
  colnames(stt_table) <- c("Time","nI","nA","nC")
  stt_table[1,] <- c(totaltime,0,0,0)

  while(timeval < totaltime)
  {
  	ext_rate <- mu * length(island_spec[,1])
  	ana_rate <- laa * length(which(island_spec[,4] == "I"))
  	clado_rate <- max(c(length(island_spec[,1]) * (lac * (1 -length(island_spec[,1])/K)),0),na.rm = T)
  	immig_rate <- max(c(mainland_n * gam * (1 - length(island_spec[,1])/K),0),na.rm = T)

  	totalrate <- ext_rate + clado_rate + ana_rate + immig_rate
  	dt <- stats::rexp(1, totalrate)

  	timeval <- timeval + dt

  	possible_event <- sample(1:4,1,replace = FALSE,c(immig_rate,ext_rate,ana_rate,clado_rate))

    ##############
    if(timeval <= totaltime)
  	{
  	  new_state <- DAISIE_sim_update_state_1_5(possible_event,maxspecID,mainland_spec,island_spec,timeval)
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
        if(!is.matrix(subset_island))
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

DAISIE_sim_update_state_1_5 <- function(possible_event,maxspecID,mainland_spec,island_spec,timeval)
{
  ##########################################
  #IMMIGRATION
  if(possible_event == 1)
  {
    colonist = DDD::sample2(mainland_spec,1)

    if(length(island_spec[,1]) != 0)
    {
      isitthere = which(island_spec[,1] == colonist)
    } else
    {
      isitthere = c()
    }

    if(length(isitthere) == 0)
    {
      island_spec = rbind(island_spec,c(colonist,colonist,timeval,"I",NA,NA,NA))
    }

    if(length(isitthere) != 0)
    {
      island_spec[isitthere,] = c(colonist,colonist,timeval,"I",NA,NA,NA)
    }
  }

  ##########################################
  #EXTINCTION
  if(possible_event == 2)
  {
    extinct = DDD::sample2(1:length(island_spec[,1]),1)
    #this chooses the row of species data to remove

    typeofspecies = island_spec[extinct,4]

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

          #first find species with same ancestor AND arrival time
          sisters = intersect(which(island_spec[,2] == island_spec[extinct,2]),which(island_spec[,3] == island_spec[extinct,3]))
          survivors = sisters[which(sisters != extinct)]

          if(length(sisters) == 2)
          {
            #survivors status becomes anagenetic
            island_spec[survivors,4] = "A"
            island_spec[survivors,c(5,6)] = c(NA,NA)
            island_spec[survivors,7] = "Clado_extinct"
            island_spec = island_spec[-extinct,]
          } else if(length(sisters) >= 3)
          {
            numberofsplits = nchar(island_spec[extinct,5])

            mostrecentspl = substring(island_spec[extinct,5],numberofsplits)

            if(mostrecentspl == "B")
            {
              sistermostrecentspl = "A"
            } else if(mostrecentspl == "A")
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
    immi_specs = which(island_spec[,4] == "I")

    #we only allow immigrants to undergo anagenesis
    if(length(immi_specs) == 1)
    {
      anagenesis = immi_specs
    } else if(length(immi_specs) > 1)
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
  return(list(island_spec = island_spec,maxspecID = maxspecID))
}
