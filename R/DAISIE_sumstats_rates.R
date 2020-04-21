#' Calculate summary statistics of per capita rates
#'
#' Calculates the summary statistics of per capita rates throught an ontogeny
#' simulation
#'
#' @inheritParams default_params_doc
#' @author Pedro Neves
#' @return A named list with mean and median values of specified parameters
#' @examples
#' pars <- c(0.01, 1, 20, 0.01, 1)
#' area_pars <- create_area_pars(
#'   max_area = 13500,
#'   current_area = 3500,
#'   proportional_peak_t = 0.1,
#'   total_island_age = 15,
#'   sea_level_amplitude = 0,
#'   sea_level_frequency = 0,
#'   island_gradient_angle = 0
#' )
#' hyper_pars <- create_hyper_pars(d = 0.2, x = 0.1)
#' island_ontogeny <- 1
#' sea_level <- 0
#' extcutoff <- 1100
#' resol <- 100
#' totaltime <- 10
#' mainland_n <- 1000
#'
#' mean_med <- DAISIE_calc_sumstats_pcrates(
#'   pars = pars,
#'   totaltime = totaltime,
#'   area_pars = area_pars,
#'   hyper_pars = hyper_pars,
#'   island_ontogeny = island_ontogeny,
#'   sea_level = sea_level,
#'   extcutoff = extcutoff,
#'   mainland_n = mainland_n,
#'   resol = resol
#' )
#' @export
DAISIE_calc_sumstats_pcrates <- function(
  pars,
  totaltime,
  area_pars,
  hyper_pars,
  island_ontogeny = 1,
  sea_level = 0,
  extcutoff = 100,
  mainland_n = 1000,
  resol = 100
) {
  peak <- calc_peak(totaltime = totaltime,
                    area_pars = area_pars)

  testit::assert(pars > 0)
  testit::assert(resol > 0)
  testit::assert(are_area_pars(area_pars))
  testit::assert(length(pars) == 5)
  testit::assert(totaltime > 0)
  testit::assert(mainland_n > 0)

  lac <- pars[1]
  mu <- pars[2]
  K <- pars[3]
  gam <- pars[4]

  # Initialize time vector given resolution and totaltime
  res <- 1 / resol
  time_vector <- seq(0, totaltime, by = res)

  A_vector <- sapply(
    X = time_vector,
    FUN = island_area,
    totaltime = totaltime,
    area_pars = area_pars,
    peak = peak,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level
  )

  # Calculate rate vectors
  clado_rates <- sapply(
    X = A_vector,
    FUN = get_clado_rate,
    lac = lac,
    hyper_pars = hyper_pars,
    num_spec = 1,
    K = K
  )

  ext_rates <- sapply(
    X = A_vector,
    FUN = get_ext_rate,
    extcutoff = extcutoff,
    mu = mu,
    hyper_pars = hyper_pars,
    num_spec = 1,
    K = K
  )

  immig_rates <- sapply(
    X = A_vector,
    FUN = get_immig_rate,
    gam = gam,
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
