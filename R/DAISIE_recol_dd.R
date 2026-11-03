#' Diversity dependence of recolonisation
#'
#' Recolonisation is the colonisation of the island by a mainland species that
#' is already present on the island as a non-endemic species. It resets the
#' colonisation time of that species but does not add a species to the
#' island. Recolonisation is therefore not
#' diversity dependent: it occurs at rate \code{gam} per mainland species
#' present as a non-endemic, whereas colonisation by a mainland species that is
#' not present as a non-endemic occurs at the (possibly diversity-dependent)
#' rate \code{gam * (1 - N / K)} (with the island area in the time-dependent
#' models). DAISIE up to version 5.0.2.1 multiplied recolonisation by the same
#' factor.
#'
#' The previous behaviour can be restored, for instance to reproduce earlier
#' results, with the argument \code{recol_dd = TRUE} of the simulation and
#' likelihood functions or for a whole session with
#' \code{options(DAISIE.recol_dd = TRUE)}. The argument sets the option for
#' the duration of the call, so functions called from within (for instance
#' the likelihood during maximum likelihood optimisation) use the same
#' setting. Without diversity dependence of colonisation (\code{ddmodel} other
#' than 11 and 21, or \code{K = Inf}) the two behaviours are identical.
#'
#' @return \code{DAISIE_recol_dd()} returns the current setting (logical).
#' @name DAISIE_recol_dd
#' @export
DAISIE_recol_dd <- function() {
  isTRUE(getOption("DAISIE.recol_dd", FALSE))
}

#' Set the diversity dependence of recolonisation for the calling function
#'
#' @param recol_dd Logical, see \code{\link{DAISIE_recol_dd}}.
#' @param envir Environment of the calling function; the previous setting is
#'   restored when it exits.
#' @return Nothing; called for its side effect.
#' @keywords internal
local_recol_dd <- function(recol_dd, envir = parent.frame()) {
  if (length(recol_dd) != 1 || is.na(recol_dd) || !is.logical(recol_dd)) {
    stop("recol_dd must be TRUE or FALSE.")
  }
  old <- options(DAISIE.recol_dd = recol_dd)
  do.call(base::on.exit, list(call("options", old), add = TRUE), envir = envir)
  invisible(NULL)
}

#' Per-capita recolonisation rate
#'
#' Rate at which a mainland species that is present on the island as a
#' non-endemic species recolonises the island.
#'
#' @param gam Colonisation rate per mainland species.
#' @param gamvec Per-capita colonisation rates including diversity dependence
#'   (as used for species that are not present as non-endemic), with the
#'   numbers of species of the states for which the recolonisation rate is
#'   needed.
#' @return \code{gamvec} if recolonisation is diversity dependent
#'   (\code{DAISIE_recol_dd()} is \code{TRUE}), otherwise \code{gam} repeated.
#' @keywords internal
get_recol_rate_per_capita <- function(gam, gamvec) {
  if (DAISIE_recol_dd()) {
    return(gamvec)
  }
  return(rep(gam, length(gamvec)))
}
