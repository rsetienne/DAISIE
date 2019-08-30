#' Create the rates per ?island
#'
#' @param immig_rate immigration rate
#' @param ext_rate extinction rate
#' @param ana_rate anagenesis rate
#' @param clado_rate cladogenesis rate
#' @param ext_rate_max maximum extinction rate
#' @param immig_rate_max maximum immigration rate
#' @param clado_rate_max maximum cladogenesis rate
#' @param Tpars A named list containing diversification rates considering two trait states:
#' \itemize{
#'   \item{[1]:A numeric with the per capita transition rate with state1}
#'   \item{[2]:A numeric with the per capita immigration rate with state2}
#'   \item{[3]:A numeric with the per capita extinction rate with state2}
#'   \item{[4]:A numeric with the per capita anagenesis rate with state2}
#'   \item{[5]:A numeric with the per capita cladogenesis rate with state2}
#'   \item{[6]:A numeric with the per capita transition rate with state2} 
#'   \item{[7]:A numeric with the number of species with trait state 2 on mainland} 
#' }
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
#'       Tpars = NULL
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
  Tpars = NULL) {
  if (is.null(Tpars)){
    rates <- list(
      immig_rate = immig_rate,
      ext_rate = ext_rate,
      ana_rate = ana_rate,
      clado_rate = clado_rate,
      ext_rate_max = ext_rate_max,
      immig_rate_max = immig_rate_max,
      clado_rate_max = clado_rate_max)
  }else{
      rates <- list(
        immig_rate = immig_rate,
        ext_rate = ext_rate,
        ana_rate = ana_rate,
        clado_rate = clado_rate,
        ext_rate_max = ext_rate_max,
        immig_rate_max = immig_rate_max,
        clado_rate_max = clado_rate_max,
        trans_rate = Tpars$trans_rate,
        immig_rate2 = Tpars$immig_rate2,
        ext_rate2 = Tpars$ext_rate2,
        ana_rate2 = Tpars$ana_rate2,
        clado_rate2 = Tpars$clado_rate2,
        trans_rate2 = Tpars$trans_rate2)
  }
  rates
}
