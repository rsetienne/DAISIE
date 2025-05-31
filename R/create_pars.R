#' Create vector of model parameters
#'
#' @inheritParams default_params_doc
#' @return Numeric vector with 5 elements
#' @export
#'
#' @examples
#' create_pars(clado_rate = 1,
#'             ext_rate = 1,
#'             carr_cap = 10,
#'             immig_rate = 0.1,
#'             ana_rate = 1)
create_pars <- function(clado_rate,
                        ext_rate,
                        carr_cap,
                        immig_rate,
                        ana_rate) {
  testit::assert(clado_rate >= 0.0)
  testit::assert(ext_rate >= 0.0)
  testit::assert(carr_cap >= 0.0)
  testit::assert(immig_rate >= 0.0)
  testit::assert(ana_rate >= 0.0)
  pars <- c(clado_rate, ext_rate, carr_cap, immig_rate, ana_rate)
  return(pars)
}

#' Create named list of area parameters
#'
#' @inheritParams default_params_doc
#'
#' @return list of numerical values containing area and sea level parameters
#' for island ontogeny simulation
#' @export
#' @author Richel J.C Bilderbeek, Joshua Lambert, Pedro Neves
#'
#' @examples
#' area_pars <- create_area_pars(
#'   max_area = 10,
#'   current_area = 1,
#'   proportional_peak_t = 0.5,
#'   total_island_age = 5,
#'   sea_level_amplitude = 5,
#'   sea_level_frequency = 10,
#'   island_gradient_angle = 0
#'
#'
#' )
create_area_pars <- function(max_area,
                             current_area,
                             proportional_peak_t,
                             total_island_age,
                             sea_level_amplitude,
                             sea_level_frequency,
                             island_gradient_angle) {
  testit::assert(max_area > 0.0)
  testit::assert(current_area > 0.0)
  testit::assert(proportional_peak_t >= 0.0)
  testit::assert(proportional_peak_t <= 1.0)
  testit::assert(total_island_age >= 0.0)
  testit::assert(sea_level_amplitude >= 0.0)
  testit::assert(sea_level_frequency >= 0.0)
  testit::assert(island_gradient_angle >= 0)
  testit::assert(island_gradient_angle <= 90)
  list(max_area = max_area,
       current_area = current_area,
       proportional_peak_t = proportional_peak_t,
       total_island_age = total_island_age,
       sea_level_amplitude = sea_level_amplitude,
       sea_level_frequency = sea_level_frequency,
       island_gradient_angle = island_gradient_angle)
}

#' Create list of hyperparameters
#'
#' @inheritParams default_params_doc
#'
#' @return Named list with hyperparameters
#' @export
#' @author Pedro Neves, Joshua Lambert
#'
#' @examples
#' hyper_pars <- create_hyper_pars(d = 0.027, x = 0.15)
create_hyper_pars <- function(d, x) {
  testit::assert(d >= 0.0)
  testit::assert(is.numeric(x))
  list(
    d = d,
    x = x
  )
}

#' Create named list of trait state parameters
#'
#' @param trans_rate   A numeric with the per capita transition rate with state1
#' @param immig_rate2  A numeric with the per capita immigration rate with state2
#' @param ext_rate2    A numeric with the per capita extinction rate with state2
#' @param ana_rate2    A numeric with the per capita anagenesis rate with state2
#' @param clado_rate2  A numeric with the per capita cladogenesis rate with state2
#' @param trans_rate2  A numeric with the per capita transition rate with state2
#' @param M2           A numeric with the number of species with trait state 2 on mainland
#'
#' @return list of numerical values containing trait state parameters
#' @export
#'
#' @examples
#' trait_pars <- create_trait_pars(
#'    trans_rate = 0.5,
#'    immig_rate2 = 0.1,
#'    ext_rate2 = 0.2,
#'    ana_rate2 = 0.3,
#'    clado_rate2 = 0.4,
#'    trans_rate2 = 0.5,
#'    M2 = 1000
#'  )
create_trait_pars <- function(trans_rate,
                              immig_rate2,
                              ext_rate2,
                              ana_rate2,
                              clado_rate2,
                              trans_rate2,
                              M2) {
  # testit::assert(is.numeric(trans_rate))
  # testit::assert(is.numeric(immig_rate2))
  # testit::assert(is.numeric(ext_rate2))
  # testit::assert(is.numeric(ana_rate2))
  # testit::assert(is.numeric(clado_rate2))
  # testit::assert(is.numeric(trans_rate2))
  # testit::assert(floor(M2) == M2)
  # testit::assert(trans_rate >= 0.0)
  # testit::assert(immig_rate2 >= 0.0)
  # testit::assert(ext_rate2 >= 0.0)
  # testit::assert(ana_rate2 >= 0.0)
  # testit::assert(clado_rate2 >= 0.0)
  # testit::assert(trans_rate2 >=0.0)
  # testit::assert(M2 >=0)
  list(trans_rate = trans_rate,
       immig_rate2 = immig_rate2,
       ext_rate2 = ext_rate2,
       ana_rate2 = ana_rate2,
       clado_rate2 = clado_rate2,
       trans_rate2 = trans_rate2,
       M2 = M2)
}



#' Create named list of trait state parameters
#'
#' @param trans_rate   A numeric with the per capita transition rate with state1
#' @param immig_rate2  A numeric with the per capita immigration rate with state2
#' @param ext_rate2    A numeric with the per capita extinction rate with state2
#' @param ana_rate2    A numeric with the per capita anagenesis rate with state2
#' @param clado_rate2  A numeric with the per capita cladogenesis rate with state2
#' @param trans_rate2  A numeric with the per capita transition rate with state2
#' @param M2           A numeric with the number of species with trait state 2 on mainland
#' @param K2           A numeric with the carrying capacity for state 2
#'
#' @return list of numerical values containing trait state parameters
#' @export
#'
create_trait_pars_2K <- function(trans_rate,
                                 immig_rate2,
                                 ext_rate2,
                                 ana_rate2,
                                 clado_rate2,
                                 trans_rate2,
                                 M2,
                                 K2) {
  list(trans_rate = trans_rate,
       immig_rate2 = immig_rate2,
       ext_rate2 = ext_rate2,
       ana_rate2 = ana_rate2,
       clado_rate2 = clado_rate2,
       trans_rate2 = trans_rate2,
       M2 = M2,
       K2 = K2)
}


#' Creates the list object for CS_version argument in DAISIE_ML_CS
#'
#' @param model the CS model to run, options are \code{1} for single rate
#' DAISIE model, \code{2} for multi-rate DAISIE, or \code{0} for IW test
#' model
#' @param function_to_optimize likelihood function that must be optimized in ML,
#'  either 'DAISIE', 'DAISIE_approx', or 'DAISIE_DE'
#' @param relaxed_par the parameter to relax (integrate over). Options are
#' \code{"cladogenesis"},
#' \code{"extinction"},
#' \code{"carrying_capacity"},
#' \code{"immigration"},
#' \code{"anagenesis"}
#' @param par_sd standard deviation of the parameter to relax
#' @param par_upper_bound upper bound of the parameter to relax
#' @param integration_method method of integration, either 'standard','stratified'
#' or 'MC'
#' @param seed seed of the random number generator in case of 'MC'
#' @param sample_size size of sample in case of 'MC' or 'stratified'
#' @param parallel use parallel computing or not in case of 'MC' or 'stratified'
#' @param n_cores number of cores to use when run in parallel
#' @return A list of four elements
#' \itemize{
#'   \item{model: the CS model to run, options are \code{1} for single rate
#'   DAISIE model, \code{2} for multi-rate DAISIE, or \code{0} for IW test
#'   model}
#'   \item{fumction_to_optimize} likelihood function that must be optimized in
#'   ML, either 'DAISIE', 'DAISIE_approx', or 'DAISIE_DE'
#'   \item{relaxed_par: the parameter to relax (integrate over), for model = 2.}
#'   \item{par_sd: the standard deviation of the parameter to relax}
#'   \item{par_upperbound: upper bound of the parameter to relax.}
#'   \item{integration_method: method of integration, either 'standard',
#'   'stratified' or 'MC'}
#'   \item{seed: random seed in case of integration_method = 'MC'}
#'   \item{sample_size: size of sample in case of integration_method = 'MC'
#'   or 'stratified'}
#'   \item{parallel: use parallel computing or not in case of integration_method
#'   = 'MC' or 'stratified'}
#'   \item{n_cores: number of cores to use when run in parallel}
#' }
#' @export
create_CS_version <- function(model = 1,
                              function_to_optimize = 'DAISIE',
                              relaxed_par = NULL,
                              par_sd = 0,
                              par_upper_bound = Inf,
                              integration_method = 'standard',
                              seed = 42,
                              sample_size = 100,
                              parallel = FALSE,
                              n_cores = 1) {

  if (model != 1 && model != 2 && model != 3) {
    stop("model must be either 1, 2 or 3")
  }
  if (model == 2 && is.null(relaxed_par)) {
    stop("relaxed_par required for multi-rate model")
  }
  if (model == 2) {
    if(integration_method == 'MC')
      CS_version <- list(model = model,
                         function_to_optimize = function_to_optimize,
                         relaxed_par = relaxed_par,
                         par_sd = par_sd,
                         par_upper_bound = par_upper_bound,
                         integration_method = integration_method,
                         seed = seed,
                         sample_size = sample_size,
                         parallel = parallel,
                         n_cores = n_cores)
    else if(integration_method == 'stratified')
      CS_version <- list(model = model,
                         function_to_optimize = function_to_optimize,
                         relaxed_par = relaxed_par,
                         par_sd = par_sd,
                         par_upper_bound = par_upper_bound,
                         integration_method = integration_method,
                         sample_size = sample_size,
                         parallel = parallel,
                         n_cores = n_cores)
    else if(integration_method == 'standard')
      CS_version <- list(model = model,
                         function_to_optimize = function_to_optimize,
                         relaxed_par = relaxed_par,
                         par_sd = par_sd,
                         par_upper_bound = par_upper_bound,
                         integration_method = 'standard')
  } else {
    CS_version <- list(model = model,
                       function_to_optimize = function_to_optimize)
  }
  return(CS_version)
}

