#' Diversity dependence of recolonization in CS_version
#'
#' Recolonization is the colonization of the island by a mainland species that
#' is already present on the island as a non-endemic species. It resets the
#' colonization time of that species but does not add a species to the island.
#' With \code{dd_recolonization = FALSE} (the default) recolonization is
#' therefore not diversity dependent: it occurs at rate \code{gam} per mainland
#' species present as a non-endemic, whereas colonization by a mainland species
#' that is not present as a non-endemic occurs at the (possibly
#' diversity-dependent) rate \code{gam * (1 - N / K)}. With
#' \code{dd_recolonization = TRUE} recolonization is multiplied by the same
#' factor, as in DAISIE up to version 5.0.3.
#'
#' @param CS_version A numeric value (the model; then recolonization is not
#'   diversity dependent) or a list created by \code{\link{create_CS_version}}.
#' @return \code{TRUE} if \code{CS_version} is a list with element
#'   \code{dd_recolonization = TRUE}, \code{FALSE} otherwise.
#' @keywords internal
get_dd_recolonization <- function(CS_version) {
  if (is.list(CS_version) && !is.null(CS_version$dd_recolonization)) {
    check_dd_recolonization(CS_version$dd_recolonization)
    return(CS_version$dd_recolonization)
  }
  return(FALSE)
}

#' Check the dd_recolonization setting
#'
#' @param dd_recolonization The value to check.
#' @return Nothing; stops with an error if \code{dd_recolonization} is not
#'   \code{TRUE} or \code{FALSE}.
#' @keywords internal
check_dd_recolonization <- function(dd_recolonization) {
  if (!is.logical(dd_recolonization) || length(dd_recolonization) != 1 ||
      is.na(dd_recolonization)) {
    stop("dd_recolonization must be TRUE or FALSE.")
  }
  invisible(NULL)
}

#' Per-capita recolonization rate
#'
#' Rate at which a mainland species that is present on the island as a
#' non-endemic species recolonizes the island.
#'
#' @param gam Colonization rate per mainland species.
#' @param gamvec Per-capita colonization rates including diversity dependence
#'   (as used for species that are not present as non-endemic), with the
#'   numbers of species of the states for which the recolonization rate is
#'   needed.
#' @param dd_recolonization Logical: is recolonization diversity dependent?
#' @return \code{gamvec} if \code{dd_recolonization} is \code{TRUE}, otherwise
#'   \code{gam} repeated.
#' @keywords internal
get_recol_rate_per_capita <- function(gam, gamvec, dd_recolonization) {
  if (dd_recolonization) {
    return(gamvec)
  }
  return(rep(gam, length(gamvec)))
}
