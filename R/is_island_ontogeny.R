#' Check if island_ontogeny is correct after user input
#'
#' @param island_ontogeny Character string that can be \code{"const"},
#' \code{"linear"} or \code{"beta"} depending on type of island ontogeny desired
#' @seealso is_island_ontogeny_runtime
#' @return Boolean stating if island_ontogeny is correct.
#' @export
is_island_ontogeny_input <- function(island_ontogeny) {
  if (class(island_ontogeny) != class(character())) return(FALSE)
  if (island_ontogeny != "const" && island_ontogeny != "linear" && island_ontogeny != "beta") return(FALSE)
  TRUE
}

#' Check if island_ontogeny is correct during runtime (i.e. numeric)
#'
#' @param island_ontogeny Character string that can be \code{"const"},
#' \code{"linear"} or \code{"beta"} depending on type of island ontogeny desired
#' @seealso is_island_ontogeny_runtime
#' @return Boolean stating if island_ontogeny is correct.
#' @export
is_island_ontogeny_runtime <- function(island_ontogeny) {
  if (class(island_ontogeny) != class(numeric())) return(FALSE)
  if (island_ontogeny != 0 && island_ontogeny != 1 && island_ontogeny != 2) return(FALSE)
  TRUE
}