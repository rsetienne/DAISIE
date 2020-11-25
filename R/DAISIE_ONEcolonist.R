#' Title
#'
#' @inheritParams default_params_doc
#'
#' @return a list with these elements:
#' \itemize{
#'   item{[1]: stt_table, the same stt_table as put in}
#'   item{[2]: branching_times, branching times}
#'   item{[3]: stac, status of colonist; see Details secion for more info}
#'   item{[4]: missing_species, ?the number of missing species}
#'   item{[5]: all_colonisations, ordered numeric vector containing all events
#'     pertaining to extant species. This includes all colonisation and
#'     branching times of extant lineages.}
#'   item{[6]: non-endemic species}
#'   item{[7]: endemic species}
#'   }
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
    all_event_times <- unique(
      sort(c(col_times, btimes_all_clado_desc), decreasing = TRUE)
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
      testit::assert(identical(all_event_times, btimes_all_clado_desc))
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

    ### create table with information on other clades with same ancestor
    ### but this information is not used in DAISIE_ML
    oldest <- which(as.numeric(island_spec[, "Colonisation time (BP)"]) ==
                      max(as.numeric(island_spec[, "Colonisation time (BP)"])))

    youngest_table <- island_spec[-oldest, ]
    if (is.character(youngest_table) && !is.matrix(youngest_table)) {
      youngest_table <- t(as.matrix(youngest_table))
    }

    uniquecol <- as.numeric(unique(youngest_table[, "Colonisation time (BP)"]))

    # all_colonisations section
    all_event_times <- c(time, all_event_times)
    descendants$all_colonisations <- all_event_times
  }
  return(descendants)
}
