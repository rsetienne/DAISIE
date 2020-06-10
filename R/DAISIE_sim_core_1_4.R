#' Internal function of the DAISIE simulation
#' Taken from CRAN, commit https://github.com/richelbilderbeek/DAISIE/commit/c700b0fcf9e2c2b5d7248f02af7596fac5a2f573#diff-ddae7ad3ae2def3cb66ecf8a8a45cc41
#' @param time simulated amount of time
#' @param mainland_n number of mainland species, that
#'   is, the number of species that can potentially colonize the island.
#'   If \code{\link{DAISIE_sim_constant_rate}()} uses a clade-specific diversity dependence,
#'   this value is set to 1.
#'   If \code{\link{DAISIE_sim_constant_rate}()} uses an island-specific diversity dependence,
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
    dt <- stats::rexp(1, totalrate)
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
      island <- DAISIE_ONEcolonist(time,island_spec,stt_table)
    }

    if(mainland_n>1) {

      ### number of colonists present
      colonists_present = sort(as.numeric(unique(island_spec[,'Mainland Ancestor'])))
      number_colonists_present = length(colonists_present)

      island_clades_info<-list()

      for (i in 1:number_colonists_present) {

        subset_island<-island_spec[which(island_spec[,'Mainland Ancestor']==colonists_present[i]),]

        if(!is.matrix(subset_island)) { subset_island<-rbind(subset_island[1:7])
        colnames(subset_island) = cnames}

        island_clades_info[[i]]<-DAISIE_ONEcolonist(time,
                                                    island_spec  =subset_island,
                                                    stt_table = NULL)
        island_clades_info[[i]]$stt_table<-NULL

      }

      island = list(stt_table = stt_table, taxon_list = island_clades_info)

    }
  }
  return(island)
}

#' Calculate the clade-wide extinction rate
#' @param ps_ext_rate per species extinction rate
#' @param n_species number of species in that clade
#' @return the clade's extinction rate
#' @author Richel J.C. Bilderbeek
#' @examples
#'   testit::assert(
#'     DAISIE_calc_clade_ext_rate(
#'       ps_ext_rate = 0.2,
#'       n_species = 4
#'     ) == 0.8
#'   )
#' @export
DAISIE_calc_clade_ext_rate <- function(ps_ext_rate, n_species) {
  testit::assert(ps_ext_rate >= 0.0)
  testit::assert(n_species >= 0)
  ps_ext_rate * n_species
}

#' Calculate the clade-wide effective anagenesis rate.
#' With 'effective', this means that if an immigrant
#' undergoes anagenesis, it will become a new species.
#' Would such a species undergo anagenesis again, no net new
#' species is created; the species only gets renamed
#' @param ps_ana_rate per species anagensis rate
#' @param n_immigrants number of immigrants in that clade
#' @return the clade's effective anagenesis rate
#' @author Richel J.C. Bilderbeek
#' @examples
#'   testit::assert(
#'     DAISIE_calc_clade_ana_rate(
#'       ps_ana_rate = 0.3,
#'       n_immigrants = 5
#'     ) == 1.5
#'   )
#' @export
DAISIE_calc_clade_ana_rate <- function(ps_ana_rate, n_immigrants) {
  testit::assert(ps_ana_rate >= 0.0)
  testit::assert(n_immigrants >= 0)
  ps_ana_rate * n_immigrants
}

#' Calculate the clade-wide cladogenesis rate.
#' @param ps_clado_rate per species cladogenesis rate
#' @param n_species number of species in that clade
#' @param carr_cap carrying capacity, number of species this clade will
#'   grow to
#' @return the clade's cladogenesis rate, which is at least zero. This
#'   rate will be zero if there are more species than the carrying capacity
#'   allows for
#' @note For clade-specific carrying capacity,
#'   each clade is simulated seperately in \code{\link{DAISIE_sim}}
#' @author Richel J.C. Bilderbeek
#' @examples
#'   testit::assert(
#'     DAISIE_calc_clade_clado_rate(
#'       ps_clado_rate = 0.2,
#'       n_species = 5,
#'       carr_cap = 10
#'     ) == 0.5
#'   )
#'   testit::assert(
#'     DAISIE_calc_clade_clado_rate(
#'       ps_clado_rate = 0.2,
#'       n_species = 2,
#'       carr_cap = 1
#'     ) == 0.0
#'   )
#' @export
DAISIE_calc_clade_clado_rate <- function(ps_clado_rate, n_species, carr_cap) {
  testit::assert(ps_clado_rate >= 0.0)
  testit::assert(n_species >= 0)
  testit::assert(carr_cap >= 0)
  return(max(
    0.0,
    n_species * ps_clado_rate * (1.0 - (n_species / carr_cap))
  ))
}

#' Calculate the clade-wide immigration rate.
#' @param ps_imm_rate per species immigration rate
#' @param n_island_species number of species in that clade on the island
#' @param n_mainland_species number of species in that clade on the mainland
#' @param carr_cap carrying capacity, number of species this clade will
#'   grow to
#' @return the clade's immigration rate, which is at least zero. This
#'   rate will be zero if there are more species than the carrying capacity
#'   allows for
#' @author Richel J.C. Bilderbeek
#' @examples
#'   testit::assert(
#'     DAISIE_calc_clade_imm_rate(
#'       ps_imm_rate = 0.1,
#'       n_island_species = 5,
#'       n_mainland_species = 2,
#'       carr_cap = 10
#'     ) == 0.1
#'   )
#'   testit::assert(
#'     DAISIE_calc_clade_imm_rate(
#'       ps_imm_rate = 0.1,
#'       n_island_species = 5,
#'       n_mainland_species = 2,
#'       carr_cap = 1
#'     ) == 0.0
#'   )
#' @export
DAISIE_calc_clade_imm_rate <- function(
  ps_imm_rate,
  n_island_species,
  n_mainland_species,
  carr_cap
) {
  testit::assert(ps_imm_rate >= 0.0)
  testit::assert(n_island_species >= 0)
  testit::assert(n_mainland_species >= 0)
  testit::assert(carr_cap >= 0)
  return(max(
    0.0,
    n_mainland_species * ps_imm_rate * (1.0 - (n_island_species / carr_cap))
  ))
}
