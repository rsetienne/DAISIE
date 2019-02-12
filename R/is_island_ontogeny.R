#' Check if island_ontogeny is correct
#'
#' @param island_ontogeny Character string that can be \code{"const"},
#' \code{"linear"} or \code{"beta"} depending on type of island ontogeny desired
#'
#' @return Boolean stating if island_ontogeny is correct.
#'
is_island_ontogeny <- function(island_ontogeny) {
  if (class(island_ontogeny) != class(character())) return(FALSE)
  if (island_ontogeny != "const" && island_ontogeny != "linear" && island_ontogeny != "beta") return(FALSE)
  TRUE
}
