#' @title Simulate islands with given parameters.
#' @description This function simulates islands with given cladogenesis, extinction, Kprime,
#' immigration and anagenesis parameters. If a single parameter set is provided
#' (5 parameters) it simulates islands where all species have the same
#' macro-evolutionary process. If two paramater sets (10 parameters) are
#' provided, it simulates islands where two different macro-evolutionary
#' processes operate, one applying to type 1 species and other to type 2
#' species.
#' 
#' Returns R list object that contains the simulated islands
#' 
#' @param time Length of the simulation in time units. For examples, if an
#' island is known to be 4 million years old, setting time = 4 will simulate 
#' the entire life span of the island; setting time = 2 will stop the simulation at
#' the mid-life of the island.
#' @param M The size of mainland pool, i.e. the number of species that can 
#' potentially colonize the island.
#' @param pars Contains the model parameters: \cr \cr \code{pars[1]}
#' corresponds to lambda^c (cladogenesis rate) \cr \code{pars[2]} corresponds
#' to mu (extinction rate) \cr \code{pars[3]} corresponds to K (clade-level
#' carrying capacity). Set K=Inf for non-diversity dependence.\cr
#' \code{pars[4]} corresponds to gamma (immigration rate) \cr \code{pars[5]}
#' corresponds to lambda^a (anagenesis rate) \cr \code{pars[6]} corresponds to
#' lambda^c (cladogenesis rate) for type 2 species \cr \code{pars[7]}
#' corresponds to mu (extinction rate) for type 2 species\cr \code{pars[8]}
#' corresponds to K (clade-level carrying capacity) for type 2 species.  Set
#' K=Inf for non-diversity dependence.\cr \code{pars[9]} corresponds to gamma
#' (immigration rate) for type 2 species\cr \code{pars[10]} corresponds to
#' lambda^a (anagenesis rate) for type 2 species\cr The elements 6:10 are
#' optional and are required only when type 2 species are included.
#' @param replicates Number of island replicates to be simulated
#' @param mainland_params mainland_params parameters for simulation mainland processes.
#'   If NULL, the mainland is assumed to be static, following the
#'   assumptions of Valente et al., 2015.
#'     Else the parameters can be created by \code{DAISIE_create_mainland_params}
#' @param divdepmodel Option divdepmodel = 'CS' runs a model with clade-specific 
#' carrying capacity, where diversity-dependence operates only within single 
#' clades, i.e. only among species originating from the same mainland colonist.
#' Option divdepmodel = 'IW' runs a model with island-wide carrying capacity,
#' where diversity-dependence operates within and among clades.
<<<<<<< HEAD
#' @param island_type Option island_type = 'oceanic' is a model equal to Valente
#' et al., 2015. island_type = 'nonoceanic' is a nonoceanic model where initial
#' species richness is non-zero determined by the nonoceanic parameters.
#' @param nonoceanic A vector of length three with: the island area as a proportion
#' of the mainland, the probability of native species being nonendemic and the 
#' size of the mainland pool.
#' @param prop_type2_pool Fractioon of mainland species that belongs to the
=======
#' @param divdep The a vector of strings to determined which parameters should
#' be diversity dependent. \code{"lac"} is cladogenesis, \code{"mu"} is extinction
#' \code{"gam"} is immigration.
#' @param prop_type2_pool Fraction of mainland species that belongs to the
>>>>>>> cc65a009840efb8907afe01a0929b4c30db839aa
#' second subset of species (type 2). Applies only when two types of species
#' are simulated (length(pars) = 10).
#' @param replicates_apply_type2 Applies only when two types of species are
#' being simulated. Default replicates_apply_type2 = TRUE runs simulations 
#' until the number of islands where a type 2 species has colonised is equal to
#' the specified number of replicates. This is recommended if prop_type2_pool
#' is small of if the rate of immigration of type two species (pars[9]) is low,
#' meaning that more replicates are needed to achieved an adequate sample size
#' of islands with type 2 species. Setting replicates_apply_type2 = FALSE
#' simulates islands up to the specified number of replicates regardless of 
#' whether type 2 species have colonised or not. 
#' @param sample_freq Specifies the number of units times should be divided by
#' for plotting purposes. Larger values will lead to plots with higher 
#' resolution, but will also run slower.
#' @param plot_sims Default = TRUE plots species-through-time (STT) plots. It
#' detects how many types of species are present. If only one type of species 
#' is present, STT is plotted for all species. If two types are present, three
#' plots are produced: STT for all, STT for type 1 and STT for type 2.
#' @param island_ontogeny a string describing the type of island ontogeny. 
#' Can be \code{"const"}, \code{"beta"} for a beta function describing area 
#' through time, or \code{"linear"} for a linear function
#' @param sea_level a string describing the type of sea level change through
#' time, can be \code{"const"}, \code{"linear_pos"}, \code{"linear_neg} for a 
#' linear positive or negative change through time respectively, or \code{"sine"}
#' for a sine wave describing the sea level oscillations.
#' @param Apars A numeric vector:
#' \itemize{
#'   \item{[1]: maximum area}
#'   \item{[2]: vale from 0 to 1 indicating where in the island's history the
#'   peak area is achieved}
#'   \item{[3]: sharpness of peak}
#'   \item{[4]: total island age}
#' }
#' @param Epars A numeric vector:
#' \itemize{
#'   \item{[1]: minimum extinction when area is at peak}
#'   \item{[2]: extinction rate when current area is 0.10 of maximum area}
#' } 
#' @param divdep The a vector of strings to determined which parameters should
#' be diversity dependent. \code{"lac"} is cladogenesis, \code{"mu"} is extinction
#' \code{"gam"} is immigration.
#' @param keep_final_state logical indicating if final state of simulation 
#' should be returned. Default is \code{FALSE}. 
#' @param stored_data output of DAISIE_sim function when run with keep_final_state.
#' If not \code{NULL}
#' @param verbose \code{Default = TRUE} Give intermediate output, also if everything
#' goes ok.
#' @param ... Any arguments to pass on to plotting functions.
#'
#' @return Each simulated dataset is an element of the list, which can be
#' called using [[x]]. For example if the object is called island_replicates,
#' the last replicates is a list in itself. The first (e.g.
#' island_replicates[[x]][[1]]) element of that list has the following 
#' components: \cr \code{$island_age} - the island age \cr Then, depending on
#' whether a distinction between types is made, we have:\cr \code{$not_present}
#' - the number of mainland lineages that are not present on the island \cr
#' or:\cr \code{$not_present_type1} - the number of mainland lineages of type 1
#' that are not present on the island \cr \code{$not_present_type2} - the
#' number of mainland lineages of type 2 that are not present on the island \cr
#' \code{$stt_all} - STT table for all species on the island (nI - number of
#' non-endemic species; nA - number of anagenetic species, nC - number of
#' cladogenetic species, present - number of independent colonisations present
#' )\cr \code{$stt_stt_type1} - STT table for type 1 species on the island -
#' only if 2 types of species were simulated (nI - number of non-endemic
#' species; nA - number of anagenetic species, nC - number of cladogenetic
#' species, present - number of independent colonisations present )\cr
#' \code{$stt_stt_type2} - STT table for type 2 species on the island - only if
#' 2 types of species were simulated (nI - number of non-endemic species; nA -
#' number of anagenetic species, nC - number of cladogenetic species, present -
#' number of independent colonisations present )\cr \code{$brts_table} - Only
#' for simulations under 'IW'. Table containing information on order of events
#' in the data, for use in maximum likelihood optimization.)\cr
#' 
#' The subsequent elements of the list each contain information on a single
#' colonist lineage on the island and has 4 components:\cr
#' \code{$branching_times} - island age and stem age of the population/species
#' in the case of Non-endemic, Non-endemic_MaxAge and Endemic anagenetic
#' species. For cladogenetic species these should be island age and branching
#' times of the radiation including the stem age of the radiation.\cr
#' \code{$stac} - the status of the colonist \cr * Non_endemic_MaxAge: 1 \cr *
#' ndemic: 2 \cr * Endemic&Non_Endemic: 3 \cr * Non_endemic: 4 \cr
#' \code{$missing_species} - number of island species that were not sampled for
#' particular clade (only applicable for endemic clades) \cr \code{$type_1or2}
#' - whether the colonist belongs to type 1 or type 2 \cr
#' @author Luis Valente and Albert Phillimore
#' @seealso \code{\link{DAISIE_format_CS}} \code{\link{DAISIE_plot_sims}}
#' @references Valente, L.M., A.B. Phillimore and R.S. Etienne (2015).
#' Equilibrium and non-equilibrium dynamics simultaneously operate in the
#' Galapagos islands. Ecology Letters 18: 844-852.
#' @keywords models
#' @examples
#'  
#' cat("
#' ## Simulate 40 islands for 4 million years, where all species have equal 
#' ## rates, and plot the species-through-time plot. Pool size 1000.
#' 
#' pars_equal = c(2.550687345,2.683454548,Inf,0.00933207,1.010073119)
#' island_replicates_equal = DAISIE_sim(
#'    time = 4,
#'    M = 1000,
#'    pars = pars_equal,
#'    replicates = 40
#'    )
#' 
#' ## Simulate 15 islands for 4 million years with two types of species (type1 
#' ## and type 2), and plot the species-through-time plot. Pool size 1000. Fraction
#' ## of type 2 species in source pool is 0.163. Function will simulate until number of islands
#' ## where type 2 species has colonised is equal to number specified in replicates.
#' 
#' pars_type1 = c(0.195442017,0.087959583,Inf,0.002247364,0.873605049)
#' pars_type2 = c(3755.202241,8.909285094,14.99999923,0.002247364,0.873605049)
#' island_replicates_2types = DAISIE_sim(
#'    time = 4,
#'    M = 1000,
#'    pars = c(pars_type1,pars_type2),
#'    replicates = 15,
#'    prop_type2_pool = 0.163
#'    )
#' ## Simulate a non-oceanic island for 4 million years, and plot the species-through-time
#' ## plot. Pool size 1000. Island area as a proportion of mainland is 0.1, proportion of
#' ## native species is 0.9.
#'  pars = c(2.550687345,2.683454548,Inf,0.00933207,1.010073119)
#'  island_replicates = DAISIE_sim(
#'    time = 4,
#'    M = 1000,
#'    pars = pars
#'    replicates = 40,
#'    island_type = 'nonoceanic'
#'    nonoceanic = c(0.1, 0.9)
#'  )
#' ")
#' 
#' @export DAISIE_sim
DAISIE_sim = function(
  time,
  M,
  pars,
  replicates,
  mainland_params = NULL,
  divdepmodel = 'CS',
<<<<<<< HEAD
  island_type = 'oceanic', # 'oceanic' = intially 0 species; 'nonoceanic' = requires nonoceanic vector
  nonoceanic = NULL,
=======
  divdep = c('lac', 'gam'), #'lac is cladogenesis, 'mu' is extinction, 'gam' is immigration,and any combination
>>>>>>> cc65a009840efb8907afe01a0929b4c30db839aa
  prop_type2_pool = NA,
  replicates_apply_type2 = TRUE,
  sample_freq = 25,
  plot_sims = TRUE,
  island_ontogeny = 'const', # const = no effect; 'linear' = linear decreasing function; 'beta' = beta function;
  sea_level = 'const', # const = no effect; 'linear_pos' = linear increasing function,
  #'linear_neg' = linear decreasing function; 'sine' = symmetric fluctuations of a sin wave
  Apars = NULL,
  Epars = NULL,
  divdep = c('lac', 'gam'), #'lac is cladogenesis, 'mu' is extinction, 'gam' is immigration,and any combination
  keep_final_state = FALSE,
  stored_data = NULL,
  verbose = TRUE,
  ...
) {
  
  testit::assert(
    "island_ontogeny is not valid input. Specify 'const',\n
    'linear' or  ' beta'", is_island_ontogeny_input(island_ontogeny)
  )
  
  #TODO: TEST island_replicates INPUT! SANITIZE STORED_DATA INPUT! ASSERT + TEST
  if (!is.null(stored_data)) {
    start_midway <- TRUE
  } else {
    start_midway <- FALSE
  }
  
  # @richelbilderbeek
  if (!is.null(mainland_params)) {
    return(
      DAISIE_sim_with_mainland(
        time = time,
        M = M,
        pars = pars,
        replicates = replicates,
        mainland_params = mainland_params,
        divdepmodel = divdepmodel,
        prop_type2_pool = prop_type2_pool,
        replicates_apply_type2 = replicates_apply_type2,
        sample_freq = sample_freq
      )
    )
  }
  
  # Classic behavior
  totaltime <- time
  island_replicates  = list()
  
  if(divdepmodel =='IW')
  {
    if(length(pars) > 5)
    {
      stop('Island-wide carrying capacity model not yet implemented for two types of mainland species')
    }
    for(rep in 1:replicates)
    {
      island_replicates[[rep]] <-DAISIE_sim_core(
        time = totaltime,
        mainland_n = M,
        pars = pars,
<<<<<<< HEAD
        island_type = island_type,
        nonoceanic = nonoceanic,
=======
        divdep = divdep,
>>>>>>> cc65a009840efb8907afe01a0929b4c30db839aa
        island_ontogeny = island_ontogeny,
        sea_level = sea_level,
        Apars = Apars,
        Epars = Epars,
        keep_final_state = keep_final_state,
        island_spec = NULL
      )
      
      if (verbose == TRUE) {
        print(paste("Island replicate ", rep, sep = ""))
      }
    }
    if(island_type == "oceanic")
    {
      island_replicates = DAISIE_format_IW_oceanic(island_replicates = island_replicates,
                                                   time = totaltime,
                                                   M = M,
                                                   sample_freq = sample_freq)
    }
    if(island_type == "nonoceanic")
    {
      island_replicates = DAISIE_format_IW_nonoceanic(island_replicates = island_replicates,
                                                      time = totaltime,
                                                      M = M,
                                                      sample_freq = sample_freq)
    }
  }
  
  if(divdepmodel == 'CS')
  {
    if(length(pars) == 5)
    {
      # Midway simulation
      
      if (!is.null(stored_data)) {
        for (rep in 1:replicates)
        {
          n_colonized_replicates <- length(stored_data[[rep]]) - 1
          colonized_island_spec <- list()
          
          for (k in 1:n_colonized_replicates) {
            colonized_island_spec[[k]] <- stored_data[[rep]][[k + 1]]$island_spec
          }
          
          island_replicates = list()
          
          # Run each clade seperately
          full_list = list()
<<<<<<< HEAD
          
          # Run midway clades
          #currently only run for oceanic DAISIE
          
          for (m_spec in 1:n_colonized_replicates)
          {
            full_list[[m_spec]] <- DAISIE_sim_core(
              time = totaltime,
              mainland_n = 1,
              pars = pars,
              island_ontogeny = island_ontogeny,
              Apars = Apars,
              Epars = Epars,
              keep_final_state = keep_final_state,
              island_spec = colonized_island_spec[[m_spec]]
            )
          }
          
          # Run empty clades that didn't get colonists
          
          for (m_spec in (n_colonized_replicates + 1):1000)
          {
            full_list[[m_spec]] <- DAISIE_sim_core(
              time = totaltime,
              mainland_n = 1,
              pars = pars,
              island_ontogeny = island_ontogeny,
              Apars = Apars,
              Epars = Epars,
              keep_final_state = keep_final_state,
              island_spec = NULL
            )
=======
          if (length(colonized_island_spec) > 0) {
            
            # Run midway clades
            for (m_spec in 1:n_colonized_replicates) 
            { 	
              full_list[[m_spec]] <- DAISIE_sim_core(
                time = totaltime,
                mainland_n = 1,
                pars = pars,
                divdep = divdep,
                island_ontogeny = island_ontogeny,
                Apars = Apars,
                Epars = Epars,
                keep_final_state = keep_final_state,
                island_spec = colonized_island_spec[[m_spec]] 
              )
            }
          } else {
            
            # Run empty clades that didn't get colonists
            for (m_spec in (n_colonized_replicates + 1):1000) 
            { 	
              full_list[[m_spec]] <- DAISIE_sim_core(
                time = totaltime,
                mainland_n = 1,
                pars = pars,
                divdep = divdep,
                island_ontogeny = island_ontogeny,
                Apars = Apars,
                Epars = Epars,
                keep_final_state = keep_final_state,
                island_spec = NULL
              )
            }
>>>>>>> cc65a009840efb8907afe01a0929b4c30db839aa
          }
          
          island_replicates[[rep]] = full_list
          
          if (verbose == TRUE) {
            print(paste("Island replicate ",rep,sep = ""))
          }
        }
        
      } else {
        
        # Simulation from empty island
        for(rep in 1:replicates)
        {
          island_replicates[[rep]] = list()
          # Run each clade seperately
          full_list = list()
          
          for(m_spec in 1:M)
          {
            full_list[[m_spec]] <- DAISIE_sim_core(
              time = totaltime,
              mainland_n = 1,
              pars = pars,
<<<<<<< HEAD
              island_type = island_type,
              nonoceanic = nonoceanic,
=======
              divdep = divdep,
>>>>>>> cc65a009840efb8907afe01a0929b4c30db839aa
              island_ontogeny = island_ontogeny,
              sea_level = sea_level,
              Apars = Apars,
              Epars = Epars,
              keep_final_state = keep_final_state,
              island_spec = NULL
            )
          }
          island_replicates[[rep]] = full_list
          if (verbose == TRUE) {
            print(paste("Island replicate ",rep,sep = ""))
          }
        }
      }
    }
    
    if(length(pars) == 10)
    {
      
      if(is.na(prop_type2_pool))
      {
        stop('prop_type2_pool (fraction of mainland species that belongs to the second subset of species) must be specified when running model with two species types')
      }
      if(island_type == "nonoceanic")
      {
        stop('nonoceanic islands cannot have two type islands')
      }
      if(replicates_apply_type2 == TRUE)
      {
        island_replicates = DAISIE_sim_min_type2(time = totaltime,
                                                 M = M,
                                                 pars = pars,
                                                 replicates = replicates,
                                                 prop_type2_pool = prop_type2_pool)
      } else
      {
        for(rep in 1:replicates)
        {
          pool2 = DDD::roundn(M * prop_type2_pool)
          pool1 = M - pool2
          lac_1 = pars[1]
          mu_1 = pars[2]
          K_1 = pars[3]
          gam_1 = pars[4]
          laa_1 = pars[5]
          lac_2 = pars[6]
          mu_2 = pars[7]
          K_2 = pars[8]
          gam_2 = pars[9]
          laa_2 = pars[10]
          full_list = list()
          
          #### species of pool1
<<<<<<< HEAD
          
          for(m_spec in 1:pool1)
          {
            full_list[[m_spec]] = DAISIE_sim_core(time = totaltime,mainland_n = 1,pars = c(lac_1,mu_1,K_1,gam_1,laa_1))
=======
          for(m_spec in 1:pool1) 
          { 	
            full_list[[m_spec]] = DAISIE_sim_core(time = totaltime,mainland_n = 1,pars = c(lac_1,mu_1,K_1,gam_1,laa_1), divdep = divdep)
>>>>>>> cc65a009840efb8907afe01a0929b4c30db839aa
            full_list[[m_spec]]$type1or2  = 1
          }
          
          #### species of pool2
          
          for(m_spec in (pool1 + 1):(pool1 + pool2))
          {
            full_list[[m_spec]] = DAISIE_sim_core(time = totaltime,
                                                  mainland_n = 1,
                                                  pars = c(lac_2,mu_2,K_2,gam_2,laa_2),
                                                  divdep = divdep)
            full_list[[m_spec]]$type1or2 = 2
          }
          
          island_replicates[[rep]] = full_list
          
          if (verbose == TRUE) {
            print(paste("Island replicate ",rep,sep = ""))
          }
        }
      }
    }
    
    if (start_midway == TRUE)
    {
      island_replicates <- DAISIE_format_CS_oceanic(
        island_replicates = island_replicates,
        time = totaltime,
        M = M,
        sample_freq = sample_freq,
        start_midway = start_midway,
        verbose = verbose
      )
    }
    
    
    if(island_type == "oceanic")
    {
      island_replicates = DAISIE_format_CS_oceanic(island_replicates = island_replicates,
                                                   time = totaltime,
                                                   M = M,
                                                   sample_freq = sample_freq)
    }
    if(island_type == "nonoceanic")
    {
      island_replicates = DAISIE_format_CS_nonoceanic(island_replicates = island_replicates,
                                                      time = totaltime,
                                                      M = M,
                                                      sample_freq = sample_freq)
    }
    
  }
  
  if(plot_sims == TRUE)
  {
    DAISIE_plot_sims(island_replicates)
  }
  return(island_replicates)
}