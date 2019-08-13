
#' Create the rates per ?island
#'
#' @param immig_rate immigration rate
#' @param ext_rate extinction rate
#' @param ana_rate anagenesis rate
#' @param clado_rate cladogenesis rate
#' @param ext_rate_max maximum extinction rate
#' @param immig_rate_max maximum immigration rate
#' @param clado_rate_max maximum cladogenesis rate
#' @param immig_rate1 immigration rate with state 1
#' @param ext_rate1 extinction rate with state 1
#' @param ana_rate1 anagenesis rate with state 1
#' @param clado_rate1 cladogenesis rate with state 1
#' @param trans_rate transition rate
#' @param trans_rate1 transition rate with state 1
#'
#' @return the rates, as can be verified by \link{are_rates}
#' @export
#'
#' @examples
#' testit::assert(
#'   DAISIE:::are_rates(
#'     create_rates(
#'       immig_rate = 0.1,
#'       ext_rate = 0.2,
#'       ana_rate = 0.3,
#'       clado_rate = 0.4,
#'       ext_rate_max = 0.5,
#'       immig_rate_max = 0.6,
#'       clado_rate_max = 0.7,
#'       immig_rate1 = 0.1,
#'       ext_rate1 = 0.2,
#'       ana_rate1 = 0.3,
#'       clado_rate1 = 0.4,
#'       trans_rate = 0.5,
#'       trans_rate1 = 0.5
#'     
#'   )
#' )
#')
create_rates <- function(
  immig_rate,
  ext_rate,
  ana_rate,
  clado_rate,
  ext_rate_max,
  immig_rate_max,
  clado_rate_max,
  immig_rate1,
  ext_rate1,
  ana_rate1,
  clado_rate1,
  trans_rate,
  trans_rate1
) {
  rates <- list(
    immig_rate = immig_rate,
    ext_rate = ext_rate,
    ana_rate = ana_rate,
    clado_rate = clado_rate,
    ext_rate_max = ext_rate_max,
    immig_rate_max = immig_rate_max,
    clado_rate_max = clado_rate_max,
    immig_rate = immig_rate,
    ext_rate = ext_rate,
    ana_rate = ana_rate,
    clado_rate = clado_rate,
    trans_rate = trans_rate,
    trans_rate1 = trans_rate1
    
    )
  rates
}