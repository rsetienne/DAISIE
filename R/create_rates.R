
#' Create the rates per ?island
#'
#' @param immig_rate immigration rate
#' @param ext_rate extinction rate
#' @param ana_rate anagenesis rate
#' @param clado_rate cladogenesis rate
#' @param ext_rate_max maximum extinction rate
#' @param immig_rate_max maximum immigration rate
#' @param clado_rate_max maximum cladogenesis rate
#'
#' @return the rates, as can be verified by \link{are_rates}
#' @export
#'
#' @examples
#'   testit::assert(
#'     are_rates(
#'       create_rates(
#'         immig_rate = 0.1,
#'         ext_rate = 0.2,
#'         ana_rate = 0.3,
#'         clado_rate = 0.4,
#'         ext_rate_max = 0.5,
#'         immig_rate_max = 0.6,
#'         clado_rate_max = 0.7
#'       )
#'     )
#'   )
create_rates <- function(
  immig_rate,
  ext_rate,
  ana_rate,
  clado_rate,
  ext_rate_max,
  immig_rate_max,
  clado_rate_max
) {
  rates <- list(
    immig_rate = immig_rate,
    ext_rate = ext_rate,
    ana_rate = ana_rate,
    clado_rate = clado_rate,
    ext_rate_max = ext_rate_max,
    immig_rate_max = immig_rate_max,
    clado_rate_max = clado_rate_max)
  rates
}