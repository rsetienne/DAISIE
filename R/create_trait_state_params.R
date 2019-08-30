#' Create named list of trait state parameters
#'
#' @param trans_rate  A numeric with the per capita transition rate with state1
#' @param immig_rate2 A numeric with the per capita immigration rate with state2
#' @param ext_rate2   A numeric with the per capita extinction rate with state2
#' @param ana_rate2   A numeric with the per capita anagenesis rate with state2
#' @param clado_rate2 A numeric with the per capita cladogenesis rate with state2
#' @param trans_rate2 A numeric with the per capita transition rate with state2
#' @param M2 A numeric with the number of species with trait state 2 on mainland
#'
#' @return list of numerical values containing trait state parameters 
#' @export
#'
#' @examples 
#' testit::assert(DAISIE:::are_trait_state_params(
#'   create_trait_state_params(
#'     trans_rate = 0.5,
#'     immig_rate2 = 0.1,
#'     ext_rate2 = 0.2,
#'     ana_rate2 = 0.3,
#'     clado_rate2 = 0.4,
#'     trans_rate2 = 0.5,
#'     M2 = 1000)) == TRUE)
create_trait_state_params <- function(trans_rate,
                                      immig_rate2,
                                      ext_rate2,
                                      ana_rate2,
                                      clado_rate2,
                                      trans_rate2,
                                      M2) {
  assertive::assert_is_a_double(trans_rate)
  assertive::assert_is_a_double(immig_rate2)
  assertive::assert_is_a_double(ext_rate2)
  assertive::assert_is_a_double(ana_rate2)
  assertive::assert_is_a_double(clado_rate2)
  assertive::assert_is_a_double(trans_rate2)
  assertive::assert_is_a_number(M2)
  testit::assert(trans_rate >= 0.0)
  testit::assert(immig_rate2 >= 0.0)
  testit::assert(ext_rate2 >= 0.0)
  testit::assert(ana_rate2 >= 0.0)
  testit::assert(clado_rate2 >= 0.0)
  testit::assert(trans_rate2 >=0.0)
  testit::assert(M2 >=0)
  list(trans_rate = trans_rate,
       immig_rate2 = immig_rate2,
       ext_rate2 = ext_rate2,
       ana_rate2 = ana_rate2,
       clado_rate2 = clado_rate2,
       trans_rate2 = trans_rate2,
       M2 = M2)
}

#' Test if list has trait state parameters
#'
#' @param trait_state_params object to be tested if conforms to trait state parameters
#'
#' @return Boolean that indicates if list conforms to expected area parameters
#' as created by \link{create_trait_state_params}
#'
#' @examples
#' testit::assert(DAISIE:::are_trait_state_params(
#'   create_trait_state_params(
#'     trans_rate = 0.5,
#'     immig_rate2 = 0.1,
#'     ext_rate2 = 0.2,
#'     ana_rate2 = 0.3,
#'     clado_rate2 = 0.4,
#'     trans_rate2 = 0.5,
#'     M2 = 1000)) == TRUE)
are_trait_state_params <- function(trait_state_params) {
  if (is.null(trait_state_params) == TRUE) return(TRUE)
  if (class(trait_state_params) != class(list())) return(FALSE)
  if (!"trans_rate" %in% names(trait_state_params)) return(FALSE)
  if (!"immig_rate2" %in% names(trait_state_params)) return(FALSE)
  if (!"ext_rate2" %in% names(trait_state_params)) return(FALSE)
  if (!"ana_rate2" %in% names(trait_state_params)) return(FALSE)
  if (!"clado_rate2" %in% names(trait_state_params)) return(FALSE)
  if (!"trans_rate2" %in% names(trait_state_params)) return(FALSE)
  if (!"M2" %in% names(trait_state_params)) return(FALSE)
  if (trait_state_params$trans_rate < 0.0) return(FALSE)
  if (trait_state_params$immig_rate2 < 0.0) return(FALSE)
  if (trait_state_params$ext_rate2 < 0.0) return(FALSE)
  if (trait_state_params$ana_rate2 < 0.0) return(FALSE)
  if (trait_state_params$clado_rate2 < 0.0) return(FALSE)
  if (trait_state_params$trans_rate2 < 0.0) return(FALSE)
  if (trait_state_params$M2 < 0.0) return(FALSE)
  TRUE
}