#' Calculate summary statistics of per capita rates
#' 
#' Calculates the summary statistics of per capita rates throught an ontogeny
#' simulation
#' @param totaltime A numeric indicating total time of simulation
#' @param resol A numeric indicating number of decimals per unit time.
#' Default is 100
#' @inherit get_ext_rate
#' @inherit DAISIE_sim_core 
#' @author Pedro Neves
#' @return A named list with mean and median values of specified parameters
#' @export
DAISIE_calc_sumstats_pcrates <- function(
  pars,
  Apars,
  Epars,
  totaltime,
  island_ontogeny = "beta",
  extcutoff = 1100,
  mainland_n = 1000,
  resol = 100
) {

  testit::assert(pars > 0)
  testit::assert(resol > 0)
  testit::assert(are_area_params(Apars))
  testit::assert(Epars[1] < Epars[2])
  testit::assert(length(pars) == 5)
  testit::assert(length(Epars) == 2)
  testit::assert(totaltime > 0)
  testit::assert(mainland_n > 0)

  lac <- pars[1]
  mu <- pars[2]
  K = pars[3]
  gam <- pars[4]

  # Initialize time vector given resolution and totaltime
  res <- 1/resol
  time_vector <- seq(0, totaltime, by = res)

  # Calculate rate vectors
  clado_rates <- sapply(
    X = time_vector,
    FUN = get_clado_rate,
    lac = lac,
    Apars = Apars,
    island_ontogeny = island_ontogeny,
    island_spec = matrix(ncol = 1),
    K = K
  )

  ext_rates <- sapply(
    X = time_vector,
    FUN = get_ext_rate,
    mu = mu, 
    extcutoff = 1100,
    Apars = Apars,
    Epars = Epars, 
    island_ontogeny = island_ontogeny,
    island_spec = matrix(ncol = 1),
    K = K
  )
  
  immig_rates <- sapply(
    X = time_vector,
    FUN = get_immig_rate,
    totaltime = totaltime,
    gam = gam,
    Apars = Apars,
    island_ontogeny = island_ontogeny,
    island_spec = matrix(ncol = 1),
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