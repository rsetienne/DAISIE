#' Call the C++ simulation engine
#'
#' @description Runs the Monte Carlo loop of one of the DAISIE simulation cores
#' in C++ (see \code{src/DAISIE_sim_engine.cpp}). The engine uses R's random
#' number generator, so \code{set.seed()} controls the result, and it draws
#' from it in the same way and order as the R implementation. Hence, the R and
#' C++ implementations of a simulation core return identical simulations for
#' the same seed. The finalisation of the island (\code{DAISIE_create_island})
#' is done in R by the calling core.
#'
#' @param cfg A list with the settings of the simulation. Element \code{mode}
#' is one of \code{"const"}, \code{"time_dep"}, \code{"trait"} or \code{"sr"}.
#'
#' @return A list with the \code{island_spec} matrix and the \code{stt_table}.
#' @noRd
run_sim_engine <- function(cfg) {
  results <- .Call("daisie_sim_engine", cfg)
  stt_table <- results$stt_table
  if (ncol(stt_table) == 4) {
    colnames(stt_table) <- c("Time", "nI", "nA", "nC")
  } else {
    colnames(stt_table) <- c("Time", "nI", "nA", "nC", "nI2", "nA2", "nC2")
  }
  list(island_spec = results$island_spec, stt_table = stt_table)
}

#' Convert the initial island state of a core to the format of the C++ engine
#'
#' @param island_spec Island specification matrix or NULL (empty island).
#' @return Character matrix with seven columns.
#' @noRd
daisie_engine_island_spec <- function(island_spec) {
  if (is.null(island_spec) || length(island_spec) == 0) {
    return(matrix(character(0), nrow = 0, ncol = 7))
  }
  if (!is.matrix(island_spec)) {
    island_spec <- rbind(island_spec)
  }
  island_spec
}
