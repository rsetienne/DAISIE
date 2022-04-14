#' DEPRECATED - Algorithm component of DAISIE_SR_sim.
#'
#'#' @details This function's use has been deprecated in favour
#' of \code{\link{DAISIE_sim_cr_shift}()}. Please use that
#' function instead.
#'
#' @param time Length of the simulation in time units. For example, if an
#' island is know to be 4 million years old, setting time = 4 will simulate
#' entire life span of the island; setting time = 2 will stop the simulation at
#' the mid-life of the island.
#' @param mainland_n A numeric stating the number of mainland species, that
#' is the number of species that can potentially colonize the island.
#' If using a clade-specific diversity dependence, this value is set to 1.
#' If using an island-wide diversity dependence, this value is set to the
#' number of mainland species.
#' @param pars Contains the model parameters: \cr \cr \code{pars[1]}
#' corresponds to lambda^c (cladogenesis rate) before the shift \cr \code{pars[2]} corresponds
#' to mu (extinction rate) before the shift \cr \code{pars[3]} corresponds to K (clade-level
#' carrying capacity) before the shift. Set K=Inf for non-diversity dependence.\cr
#' \code{pars[4]} corresponds to gamma (immigration rate) before the shift \cr \code{pars[5]}
#' corresponds to lambda^a (anagenesis rate) before the shift \cr \code{pars[6]} corresponds to
#' lambda^c (cladogenesis rate) after the shift \cr \code{pars[7]}
#' corresponds to mu (extinction rate) after the shift \cr \code{pars[8]}
#' corresponds to K (clade-level carrying capacity) after the shift.  Set
#' K=Inf for non-diversity dependence. \cr \code{pars[9]} corresponds to gamma
#' (immigration rate)  after the shift \cr \code{pars[10]} corresponds to
#' lambda^a (anagenesis rate) after the shift \cr \code{pars[11]} corresponds to
#' the time of shift. This is defined as time before the end of the simulation. For example,
#' setting time = 4 and pars[11] = 1.5 will simulate with pars[1:5] from 4 to 1.5 and
#' with pars[6:10] from 1.5 to 0.
#' @author Luis Valente, Albert Phillimore, and Torsten Hauffe
#' @references Hauffe, T., D. Delicado, R.S. Etienne and L. Valente (2020).
#' Lake expansion increases equilibrium diversity via the target effect of
#' island biogeography
#'
#' @return List with DAISIE simulation.
#' @keywords internal
DAISIE_SR_sim_core <- function(time,mainland_n,pars)
{
  total_time <- time
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
  stt_table[1,] <- c(total_time,0,0,0)

  while(timeval < total_time)
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
    if(timeval <= total_time)
    {
      new_state <- DAISIE_sim_update_state_cr(timeval = timeval,
                                              total_time = total_time,
                                              possible_event = possible_event,
                                              maxspecID = maxspecID,
                                              mainland_spec = mainland_spec,
                                              island_spec = island_spec,
                                              stt_table = stt_table)
      island_spec <- new_state$island_spec
      maxspecID <- new_state$maxspecID
    }
    stt_table <- rbind(stt_table,
                       c(total_time - timeval,
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
    island <- list(stt_table = stt_table, branching_times = total_time, stac = 0, missing_species = 0)
  } else
  {
    cnames <- c("Species","Mainland Ancestor","Colonisation time (BP)",
                "Species type","branch_code","branching time (BP)","Anagenetic_origin")
    colnames(island_spec) <- cnames

    ### set ages as counting backwards from present
    island_spec[,"branching time (BP)"] <- total_time - as.numeric(island_spec[,"branching time (BP)"])
    island_spec[,"Colonisation time (BP)"] <- total_time - as.numeric(island_spec[,"Colonisation time (BP)"])

    if(mainland_n == 1)
    {
      island <- DAISIE_ONEcolonist(total_time,island_spec,stt_table)
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
        island_clades_info[[i]] <- DAISIE_ONEcolonist(total_time,island_spec=subset_island,stt_table=NULL)
        island_clades_info[[i]]$stt_table <- NULL
      }
      island <- list(stt_table = stt_table, taxon_list = island_clades_info)
    }
  }
  return(island)
}
