#' Calculate summary statistics of per capita rates
#'
#' Calculates the summary statistics of per capita rates throught an ontogeny
#' simulation
#'
#' @param totaltime A numeric indicating total time of simulation.
#' @param resol A numeric > 0 indicating number of decimals per unit time.
#' @param pars A numeric vector:
#' \itemize{
#'   \item{[1]: cladogenesis rate}
#'   \item{[2]: extinction rate}
#'   \item{[3]: carrying capacity}
#'   \item{[4]: immigration rate}
#'   \item{[5]: anagenesis rate}
#' }
#' @param area_pars a named list containing area and sea level parameters as
#' created by \code{\link{create_area_pars}}:
#' \itemize{
#'   \item{[1]: maximum area}
#'   \item{[2]: value from 0 to 1 indicating where in the island's history the
#'   peak area is achieved}
#'   \item{[3]: sharpness of peak}
#'   \item{[4]: total island age}
#'   \item{[5]: amplitude of area fluctuation from sea level}
#'   \item{[6]: frequency of sine wave of area change from sea level}
#' }
#' @param ext_pars A numeric vector:
#' \itemize{
#'   \item{[1]: minimum extinction when area is at peak}
#'   \item{[2]: extinction rate when current area is 0.10 of maximum area}
#' }
#' @param island_ontogeny a numeric describing the type of island ontogeny.
#' Can be \code{0} for constant, \code{1} for a beta function describing area.
#' @param sea_level a numeric describing sea level.
#' 0 corresponds to no sea level dynamics.
#' @param extcutoff A numeric with the cutoff for extinction rate
#' preventing it from being too
#' large and slowing down simulation. Should be big. Default is 1100.
#' @param mainland_n A numeric stating the number of mainland species, that
#' is the number of species that can potentially colonize the island.
#' If DAISIE_sim function uses a clade-specific diversity dependence,
#' this value is set to 1.
#' If DAISIE_sim function uses an island-wide diversity dependence,
#' this value is set to the number of mainland species. Default is 1000
#' @inherit get_ext_rate
#' @inherit DAISIE_sim_core_time_dependent
#' @author Pedro Neves
#' @return A named list with mean and median values of specified parameters
#' @examples
#' pars <- c(0.01, 1, 20, 0.01, 1)
#' ext_pars <- c(0.2, 10)
#' area_pars <- create_area_pars(
#'   max_area = 13500,
#'   proportional_peak_t = 0.1,
#'   peak_sharpness = 1,
#'   total_island_age = 15,
#'   sea_level_amplitude = 0,
#'   sea_level_frequency = 0,
#'   island_gradient_angle = 0
#' )
#' hyper_pars <- NULL
#' dist_pars <- create_dist_pars(1)
#' island_ontogeny <- 1
#' sea_level <- 0
#' extcutoff <- 1100
#' resol <- 100
#' totaltime <- 10
#' mainland_n <- 1000
#'
#' mean_med <- DAISIE_calc_sumstats_pcrates(
#'   pars = pars,
#'   ext_pars = ext_pars,
#'   totaltime = totaltime,
#'   area_pars = area_pars,
#'   hyper_pars = hyper_pars,
#'   dist_pars = dist_pars,
#'   island_ontogeny = island_ontogeny,
#'   sea_level = sea_level,
#'   extcutoff = extcutoff,
#'   mainland_n = mainland_n,
#'   resol = resol
#' )
#' @export
DAISIE_calc_sumstats_pcrates <- function(
  pars,
  ext_pars,
  totaltime,
  area_pars = NULL,
  hyper_pars = NULL,
  dist_pars = NULL,
  island_ontogeny = 1,
  sea_level = 0,
  extcutoff = 100,
  mainland_n = 1000,
  resol = 100
) {

  meta_pars <- create_default_pars(island_ontogeny = island_ontogeny,
                                sea_level = sea_level,
                                area_pars = area_pars,
                                hyper_pars = hyper_pars,
                                dist_pars = dist_pars,
                                ext_pars = ext_pars,
                                totaltime = totaltime,
                                pars = pars)
  area_pars <- meta_pars$area_pars
  hyper_pars = meta_pars$hyper_pars
  dist_pars = meta_pars$dist_pars
  ext_pars = meta_pars$ext_pars

  testit::assert(pars > 0)
  testit::assert(resol > 0)
  testit::assert(are_area_pars(area_pars))
  testit::assert(ext_pars[1] < ext_pars[2])
  testit::assert(length(pars) == 5)
  testit::assert(length(ext_pars) == 2)
  testit::assert(totaltime > 0)
  testit::assert(mainland_n > 0)

  lac <- pars[1]
  mu <- pars[2]
  K <- pars[3]
  gam <- pars[4]

  # Initialize time vector given resolution and totaltime
  res <- 1 / resol
  time_vector <- seq(0, totaltime, by = res)

  # Calculate rate vectors
  clado_rates <- sapply(
    X = time_vector,
    FUN = get_clado_rate,
    lac = lac,
    hyper_pars = hyper_pars,
    area_pars = area_pars,
    dist_pars = dist_pars,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    num_spec = 1,
    K = K
  )

  ext_rates <- sapply(
    X = time_vector,
    FUN = get_ext_rate,
    extcutoff = extcutoff,
    hyper_pars = hyper_pars,
    area_pars = area_pars,
    ext_pars = ext_pars,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    num_spec = 1,
    K = K
  )

  immig_rates <- sapply(
    X = time_vector,
    FUN = get_immig_rate,
    totaltime = totaltime,
    gam = gam,
    hyper_pars = hyper_pars,
    area_pars = area_pars,
    dist_pars = dist_pars,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    num_spec = 1,
    mainland_n = mainland_n,
    K = K
  )

  # Calculate summary statistics
  mean_lambda_c <- mean(clado_rates)
  med_lambda_c <- stats::median(clado_rates)
  mean_mu <- mean(ext_rates)
  med_mu <- stats::median(ext_rates)
  mean_gamma <- mean(immig_rates)
  med_gamma <- stats::median(immig_rates)


  # Store in named list and return
  out <- list(
    medians = c(
    med_lambda_c = med_lambda_c,
    med_mu = med_mu,
    med_gamma = med_gamma),
    means = c(
    mean_lambda_c = mean_lambda_c,
    mean_mu = mean_mu,
    mean_gamma = mean_gamma
    )
  )
  testit::assert(is_numeric_list(out))
  testit::assert(length(out) == 2)
  testit::assert(length(out$medians) == 3)
  testit::assert(length(out$means) == 3)
  return(out)
}
