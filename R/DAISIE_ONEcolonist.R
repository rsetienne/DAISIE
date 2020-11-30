#' Convert intermediate output to final simulation output
#'
#' @inheritParams default_params_doc
#'
#' @return a list with these elements:
#' \itemize{
#'   \item{[1]: \code{stt_table}, the same stt_table as put in.}
#'   \item{[2]: \code{branching_times}, a sorted numeric vector, as required
#'     by the ML estimation functions. The first element always refers to
#'     the island age. Subsequent elements refer to colonisation, speciation and
#'     recolonisation times. The most recent recolonisation time, if any is
#'     always omitted to approximate simulation results to the mathematical
#'     formulation of the likelihood functions used for MLE.}
#'   \item{[3]: \code{stac}, status of colonist. In this function it can be
#'     returned as either 2, 4 or 3. If \code{stac} is 2, then there is only one
#'     independent colonisation present on the island and the extant species are
#'     endemic. If stac is 4, then only a singleton endemic is present at the
#'     present. If stac is 3, then recolonisation occurred, and more than one
#'     colonising lineage.}
#'   \item{[4]: \code{missing_species}, a numeric value with the number of
#'     missing species, that is, species not sampled in the phylogeny but
#'     present on the island. As this code only runs for simulation models,
#'     here \code{missing_species} is always set to 0.}
#'   \item{[5]:
#'   \code{all_colonisations}, on recolonising lineages only. It is comprised of
#'     \code{$event_times} and \code{$species_type}:
#'     \describe{
#'       \item{\code{$event_times}}{ordered numeric vectors containing all
#'       events for each extant recolonising lineage. This includes all
#'       colonisation and branching times. Each vector pertains to one
#'       colonising lineage.}
#'       \item{\code{$species_type}}{a string. Can be \code{"A"}, \code{"C"} or
#'       \code{"I"} depending on whether the extant clade is of anagenetic,
#'       cladogenetic or immigrant origin, respectively.}
#'     }
#'   }
#' }
#' @keywords internal
DAISIE_ONEcolonist <- function(time,
                               island_spec,
                               stt_table) {
  ### number of independent colonisations
  uniquecolonisation <- as.numeric(unique(
    island_spec[, "Colonisation time (BP)"]))
  number_colonisations <- length(uniquecolonisation)
  ### if there is only one independent colonisation - anagenetic and
  ### cladogenetic species are classed as stac=2; immigrant classed as stac=4:
  if (number_colonisations == 1) {
    if (island_spec[1, "Species type"] == "I") {
      descendants <- list(stt_table = stt_table,
                          branching_times = c(
                            time,
                            as.numeric(island_spec[1, "Colonisation time (BP)"])
                          ),
                          stac = 4,
                          missing_species = 0)
    }
    if (island_spec[1, "Species type"] == "A") {
      descendants <- list(stt_table = stt_table,
                          branching_times = c(
                            time,
                            as.numeric(island_spec[1, "Colonisation time (BP)"])
                          ),
                          stac = 2,
                          missing_species = 0)
    }
    if (island_spec[1, "Species type"] == "C") {
      descendants <- list(stt_table = stt_table,
                          branching_times = c(
                            time,
                            sort(
                              as.numeric(island_spec[, "branching time (BP)"]),
                              decreasing = TRUE
                            )
                          ),
                          stac = 2,
                          missing_species = 0)
    }
  }

  ### if there are two or more independent colonisations, all species are
  ### classed as stac=3 and put within same list item:
  if (number_colonisations > 1) {
    descendants <- list(stt_table = stt_table,
                        branching_times = NA,
                        stac = 3,
                        missing_species = 0,
                        all_colonisations = list())

    # Get branching and colonisation times
    btimes_all_clado_desc <- rev(
      sort(as.numeric(island_spec[, "branching time (BP)"]))
    )
    col_times <- sort(
      unique(as.numeric(island_spec[, "Colonisation time (BP)"])),
      decreasing = TRUE
    )

    # If there are endemic descendants find youngest col time
    if (length(btimes_all_clado_desc) != 0) {
      # Ensure all col_times are in b_times at this point.
      # Covers cases of one recolonization followed by cladogenesis and
      # potential extinction
      if (any(!(col_times %in% btimes_all_clado_desc))) {
        miss_col_time <- which(!(col_times %in% btimes_all_clado_desc))
        btimes_all_clado_desc <- sort(
          c(btimes_all_clado_desc, col_times[miss_col_time]),
          decreasing = TRUE
        )
      }
      youngest_col_time <- min(col_times)
      i_youngest_col_btimes <- which(btimes_all_clado_desc == youngest_col_time)

      # Remove youngest col time in branching times
      testit::assert(youngest_col_time %in% btimes_all_clado_desc)
      btimes_all_clado_desc <- btimes_all_clado_desc[-i_youngest_col_btimes]

      descendants$branching_times <- c(time, btimes_all_clado_desc)
      testit::assert(!(youngest_col_time %in% btimes_all_clado_desc))

      # If no cladogenetic species is present, remove the youngest col time
    } else if (length(btimes_all_clado_desc) == 0) {
      youngest_col_time <- min(col_times)
      i_youngest_col_time <- which(col_times == youngest_col_time)
      col_times <- col_times[-i_youngest_col_time]

      descendants$branching_times <- c(time, col_times)
    }


    # all_colonisations section
    uniquecol <- sort(as.numeric(
      unique(island_spec[, "Colonisation time (BP)"])), decreasing = TRUE
    )
    for (i in seq_along(uniquecol)) {
      descendants$all_colonisations[[i]] <- list(
        event_times = NA,
        species_type = NA
      )

      samecolonisation <- which(as.numeric(
        island_spec[, "Colonisation time (BP)"]) == uniquecol[i]
      )

      if (island_spec[samecolonisation[1], "Species type"] == "I") {
        descendants$all_colonisations[[i]]$event_times <- as.numeric(
          c(time,island_spec[samecolonisation, "Colonisation time (BP)"])
        )
        descendants$all_colonisations[[i]]$species_type <- "I"
      }

      if (island_spec[samecolonisation[1], "Species type"] == "A") {
        descendants$all_colonisations[[i]]$event_times <- as.numeric(
          c(time, island_spec[samecolonisation, "Colonisation time (BP)"])
        )
        descendants$all_colonisations[[i]]$species_type <- "A"
      }

      if (island_spec[samecolonisation[1], "Species type"] == "C") {
        descendants$all_colonisations[[i]]$event_times <-
          sort(c(time, as.numeric(
            island_spec[samecolonisation, "branching time (BP)"]
          )), decreasing = TRUE)
        descendants$all_colonisations[[i]]$species_type <- "C"
      }
    }
  }
  return(descendants)
}
