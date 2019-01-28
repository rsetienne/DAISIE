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
DAISIE_sim_core_checked <- function(
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
  DAISIE_sim_core(
    time = sim_time,
    mainland_n = n_mainland_species,
    pars = c(clado_rate, ext_rate, carr_cap, imm_rate, ana_rate)
  )
}
