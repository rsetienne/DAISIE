#' Internal function of the DAISIE simulation
#' Taken from CRAN, commit https://github.com/richelbilderbeek/DAISIE/commit/c700b0fcf9e2c2b5d7248f02af7596fac5a2f573#diff-ddae7ad3ae2def3cb66ecf8a8a45cc41
#' @param time simulated amount of time
#' @param mainland_n number of mainland species, that
#'   is, the number of species that can potentially colonize the island.
#'   If \code{\link{DAISIE_sim}} uses a clade-specific diversity dependence,
#'   this value is set to 1. 
#'   If \code{\link{DAISIE_sim}} uses an island-specific diversity dependence,
#'   this value is set to the number of mainland species.
#' @param pars a numeric vector:
#' \itemize{
#'   \item{[1]: cladogenesis rate}
#'   \item{[2]: extinction rate}
#'   \item{[3]: carrying capacity}
#'   \item{[4]: immigration rate}
#'   \item{[5]: anagenesis rate}
#' }
DAISIE_sim_core_1_4 = function(time, mainland_n, pars)
{
  lac = pars[1]
  mu = pars[2]
  K = pars[3]
  gam = pars[4]
  laa = pars[5]
  
  if(pars[4]==0) 
  {
    stop('Rate of colonisation is zero. Island cannot be colonised.')
  }  
  
  timeval = 0
  
  mainland_spec = seq(1,mainland_n,1)
  maxspecID = mainland_n
  
  island_spec = c()
  stt_table = matrix(ncol=4)
  colnames(stt_table) = c("Time","nI","nA","nC")
  stt_table[1,] = c(time,0,0,0)
  
  while(timeval < time)
  {
    n_island_species <- length(island_spec[,1])
    n_immigrants <- length(which(island_spec[,4] == "I"))
    ext_rate <- DAISIE_calc_clade_ext_rate(ps_ext_rate = mu, n_species = n_island_species)
    ana_rate <- DAISIE_calc_clade_ana_rate(ps_ana_rate = laa, n_immigrants = n_immigrants)
    clado_rate <- DAISIE_calc_clade_clado_rate(ps_clado_rate = lac, n_species = n_island_species, carr_cap = K)
    immig_rate <- DAISIE_calc_clade_imm_rate(ps_imm_rate = gam, n_island_species = n_island_species, n_mainland_species = mainland_n, carr_cap = K)
    totalrate <- ext_rate + clado_rate + ana_rate + immig_rate
    dt <- rexp(1, totalrate)
    timeval <- timeval  + dt
    possible_event <- sample(1:4,1,replace=FALSE,c(immig_rate,ext_rate,ana_rate,clado_rate))
    ##############
    if(timeval <= time)
    {  
      ##########################################
      #IMMIGRATION
      if(possible_event == 1)
      {  	
        colonist = DDD::sample2(mainland_spec,1)
        
        if(length(island_spec[,1]) != 0){isitthere = which(island_spec[,1] == colonist)}
        
        if(length(island_spec[,1]) == 0) {isitthere = c()}
        
        if(length(isitthere) == 0){island_spec = rbind(island_spec,c(colonist,colonist,timeval,"I",NA,NA,NA))}
        
        if(length(isitthere) != 0){ island_spec[isitthere,] = c(colonist,colonist,timeval,"I",NA,NA,NA)}
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
    }
    stt_table = rbind(stt_table,c(time - timeval,length(which(island_spec[,4] == "I")),
                                  length(which(island_spec[,4] == "A")),length(which(island_spec[,4] == "C"))))
  }
  
  stt_table[nrow(stt_table),1] = 0
  ############# 
  ### if there are no species on the island branching_times = island_age, stac = 0, missing_species = 0 
  if(length(island_spec[,1])==0)
  {
    island = list(stt_table = stt_table, branching_times = time, stac = 0, missing_species = 0) }
  else{
    
    cnames <- c("Species","Mainland Ancestor","Colonisation time (BP)",
                "Species type","branch_code","branching time (BP)","Anagenetic_origin")
    colnames(island_spec) <- cnames
    
    ### set ages as counting backwards from present
    island_spec[,"branching time (BP)"] = time - as.numeric(island_spec[,"branching time (BP)"])
    island_spec[,"Colonisation time (BP)"] = time - as.numeric(island_spec[,"Colonisation time (BP)"])
    
    if(mainland_n==1) {
      island <- DAISIE_ONEcolonist(time,island_spec,stt_table, keep_final_state = FALSE)
    }
    
    if(mainland_n>1) {  
      
      ### number of colonists present
      colonists_present = sort(as.numeric(unique(island_spec[,'Mainland Ancestor'])))
      number_colonists_present = length(colonists_present) 
      
      island_clades_info<-list()  
      
      for (i in 1:number_colonists_present) {
        
        subset_island<-island_spec[which(island_spec[,'Mainland Ancestor']==colonists_present[i]),] 
        
        if(class(subset_island)!='matrix') { subset_island<-rbind(subset_island[1:7])
        colnames(subset_island) = cnames}
        
        island_clades_info[[i]]<-DAISIE_ONEcolonist(time,
                                                    island_spec=subset_island,
                                                    stt_table=NULL, 
                                                    keep_final_state = FALSE)
        island_clades_info[[i]]$stt_table<-NULL
        
      }
      
      island = list(stt_table = stt_table, taxon_list = island_clades_info)
      
    }
  }
  return(island) 
}


#' Does something
#'
#' @param time simulated amount of time
#' @param island_spec matrix with current state of simulation
#' @param keep_final_state logical indicating if final state of simulation 
#' should be returned. Default is \code{FALSE}
#' @param stt_table ?Species-Through-Time table
#'
#' @return a list with these elements:
#' \itemize{
#'   item{[1]: stt_table, the same stt_table as put in}
#'   item{[2]: branching_times, branching times}
#'   item{[3]: stac, ?statuses}
#'   item{[4]: missing_species, ?the number of missing species}
#'   item{[5]: other_clades_same_ancestor, ?no idea}
#' }
DAISIE_ONEcolonist = function(time,island_spec,stt_table, keep_final_state = FALSE)
{
  
  ### number of independent colonisations
  uniquecolonisation = as.numeric(unique(island_spec[,"Colonisation time (BP)"]))
  number_colonisations = length(uniquecolonisation)
  
  ### if there is only one independent colonisation - anagenetic and cladogenetic
  #species are classed as stac=2; immigrant classed as stac=4:
  if(number_colonisations == 1)
  {
    if (island_spec[1,"Species type"] == "I")
    {
      descendants = list(stt_table = stt_table, branching_times = c(time,as.numeric(island_spec[1,"Colonisation time (BP)"])),
                         stac = 4, missing_species = 0)
    }
    if (island_spec[1,"Species type"] == "A")
    {
      descendants = list(stt_table = stt_table, branching_times = c(time,as.numeric(island_spec[1,"Colonisation time (BP)"])),
                         stac = 2,missing_species = 0)
    }
    if (island_spec[1,"Species type"] == "C")
    {
      descendants = list(stt_table = stt_table, branching_times = c(time,rev(sort(as.numeric(island_spec[,"branching time (BP)"])))),
                         stac = 2,missing_species = 0)
    }
  }
  
  ### if there are two or more independent colonisations, all species are classed as stac=3 and put within same list item:
  if(number_colonisations > 1)
  {
    descendants = list(stt_table = stt_table, branching_times = NA,stac = 3,missing_species = 0,
                       other_clades_same_ancestor = list())
    
    btimes_all_clado_desc = rev(sort(as.numeric(island_spec[,'branching time (BP)'])))
    
    if(length(btimes_all_clado_desc)!=0) { descendants$branching_times= c(time, btimes_all_clado_desc)}
    if(length(btimes_all_clado_desc)==0) { descendants$branching_times= c(time, max(as.numeric(island_spec[,"Colonisation time (BP)"])))}
    
    ### create table with information on other clades with same ancestor, but this information is not used in DAISIE_ML
    oldest = which(as.numeric(island_spec[,"Colonisation time (BP)"]) == max(as.numeric(island_spec[,"Colonisation time (BP)"])))
    
    youngest_table = island_spec[-oldest,]
    if(class(youngest_table)=='character')
    {
      youngest_table = t(as.matrix(youngest_table))
    }
    
    uniquecol = as.numeric(unique(youngest_table[,"Colonisation time (BP)"]))
    
    for(colonisation in 1:length(uniquecol))
    {
      descendants$other_clades_same_ancestor[[colonisation]] = list(brts_miss = NA,species_type = NA)
      
      samecolonisation = which(as.numeric(youngest_table[,"Colonisation time (BP)"]) == uniquecol[colonisation])
      
      if (youngest_table[samecolonisation[1],"Species type"] == "I")
      {
        descendants$other_clades_same_ancestor[[colonisation]]$brts_miss = as.numeric(youngest_table[samecolonisation,"Colonisation time (BP)"])
        descendants$other_clades_same_ancestor[[colonisation]]$species_type = "I"
      }
      
      if (youngest_table[samecolonisation[1],"Species type"] == "A")
      {
        descendants$other_clades_same_ancestor[[colonisation]]$brts_miss =  as.numeric(youngest_table[samecolonisation,"Colonisation time (BP)"])
        descendants$other_clades_same_ancestor[[colonisation]]$species_type = "A"
      }
      
      if (youngest_table[samecolonisation[1],"Species type"] == "C")
      {
        descendants$other_clades_same_ancestor[[colonisation]]$brts_miss = rev(sort(as.numeric(youngest_table[samecolonisation,"branching time (BP)"])))
        descendants$other_clades_same_ancestor[[colonisation]]$species_type = "C"
      }
    }
  }
  #### ADDS island_spec ####
  if (keep_final_state == TRUE) {
    descendants$island_spec <- island_spec
  }
  return(descendants)
}

#' Runs one DAISIE simulation with a clade-specific carrying capacity.
#' Version of \code{DAISIE_sim_core} that checks all its inputs
#' and uses descriptively named arguments
#' @param sim_time length of the simulated time
#' @param n_mainland_species number of mainland species
#' @param clado_rate cladogenesis rate 
#' @param ext_rate extinction rate
#' @param carr_cap carrying capacity
#' @param imm_rate immigration rate
#' @param ana_rate anagenesis rate
#' @return a list with these elements:
#' \describe{
#'   \item{stt_table}{a species-through-time table}
#'   \item{branching_times}{branching times}
#'   \item{stac}{
#'     the status of the colonist
#'     \itemize{
#'       \item{1: \code{Non_endemic_MaxAge} (?immigrant is present but has not formed an extant clade)}
#'       \item{2: \code{Endemic} (?immigrant is not present but has formed an extant clade)}
#'       \item{3: \code{Endemic&Non_Endemic} (?immigrant is present and has formed an extant clade)}
#'       \item{4: \code{Non_endemic} (?immigrant is present but has not formed an extant clade, and it is known when it immigrated)}
#'     }
#'   }
#'   \item{missing_species}{number of missing species}
#'   \item{other_clades_same_ancestor}{(not always present) ?no idea}
#' }
#' @author Richel J.C. Bilderbeek
DAISIE_sim_core_checked_1_4 <- function(
  sim_time, 
  n_mainland_species, 
  clado_rate, 
  ext_rate,
  carr_cap,
  imm_rate,
  ana_rate
) {
  testit::assert(sim_time > 0.0)
  testit::assert(n_mainland_species > 0)
  testit::assert(clado_rate >= 0.0)
  testit::assert(ext_rate >= 0.0)
  testit::assert(carr_cap > 0)
  testit::assert(imm_rate > 0.0)
  testit::assert(ana_rate >= 0.0)
  DAISIE_sim_core_1_4(
    time = sim_time,
    mainland_n = n_mainland_species,
    pars = c(clado_rate, ext_rate, carr_cap, imm_rate, ana_rate)
  )
}
