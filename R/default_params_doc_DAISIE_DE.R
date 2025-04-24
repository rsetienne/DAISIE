#' @title Default Parameters Documentation
#' @description This file defines the default parameters for DAISIE_DE functions.
#' @name default_params_doc
NULL

#' @param datalist A list containing colonization and branching information for island lineages.
#' This object can be created using the \code{DAISIE_dataprep()} function or manually constructed.
#' It should be a list with the following structure:
#' \itemize{
#'   \item \code{datalist[[1]]$island_age}: Age of the island.
#'   \item \code{datalist[[1]]$not_present} or (for trait-dependent cases)
#'         \code{datalist[[1]]$not_present_type1} and \code{datalist[[1]]$not_present_type2}:
#'         Number of mainland species not present on the island.
#'   \item Each subsequent element of the list corresponds to a single colonist lineage and includes:
#'     \itemize{
#'       \item \code{$colonist_name}: Name of the species or clade.
#'       \item \code{$branching_times}: A numeric vector starting with the island age, followed by colonization and speciation times.
#'       \item \code{$stac}: Colonist status, one of the following:
#'         \enumerate{
#'           \item Non_endemic_MaxAge: 1
#'           \item Endemic: 2
#'           \item Endemic & Non_Endemic: 3
#'           \item Non_Endemic: 4
#'           \item Endemic_Singleton_MaxAge: 5
#'           \item Endemic_Clade_MaxAge: 6
#'           \item Endemic & Non_Endemic_Clade_MaxAge: 7
#'         }
#'       \item \code{$missing_species}: Number of missing species for the clade (applies to endemic clades only).
#'       \item \code{$type1or2}: Lineage type (1 or 2), used in trait-dependent models.
#'     }
#' }
#'
#' @param brts The branching times of the lineage being considered in the dataset.
#' @param missnumspec The number of missing species in the lineage being considered.
#' @param pars1 A numeric vector of model parameters:
#' \itemize{
#'   \item \code{pars1[1]}: \eqn{\lambda^c} (Cladogenesis rate)
#'   \item \code{pars1[2]}: \eqn{\mu_E} (Extinction rate of endemic lineages)
#'   \item \code{pars1[3]}: \eqn{\mu_{NE}} (Extinction rate of non-endemic lineages)
#'   \item \code{pars1[4]}: \eqn{\gamma} (Colonization rate)
#'   \item \code{pars1[5]}: \eqn{\lambda^a} (Anagenesis rate)
#' }
#'
#' @param methode The numerical method to use for solving the system of differential equations.
#' @param reltolint Relative tolerance for numerical integration.
#' @param abstolint Absolute tolerance for numerical integration.
#'
