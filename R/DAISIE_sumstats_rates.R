#' Calculate summary statistics of per capita rates
#'
#' Calculates the summary statistics of per capita rates throught an ontogeny
#' simulation
#'
#' @inheritParams default_params_doc
#' @author Pedro Neves
#' @noRd
#' @return A named list with mean and median values of specified parameters
DAISIE_calc_sumstats_pcrates <- function(
  pars,
  total_time,
  area_pars,
  hyper_pars,
  island_ontogeny = 1,
  sea_level = 0,
  extcutoff = 100,
  mainland_n = 1000,
  resol = 100
) {
  peak <- calc_peak(total_time = total_time,
                    area_pars = area_pars)

  testit::assert(pars > 0)
  testit::assert(resol > 0)
  testit::assert(are_area_pars(area_pars))
  testit::assert(length(pars) == 5)
  testit::assert(total_time > 0)
  testit::assert(mainland_n > 0)

  lac <- pars[1]
  mu <- pars[2]
  K <- pars[3]
  gam <- pars[4]

  # Initialize time vector given resolution and total_time
  res <- 1 / resol
  time_vector <- seq(0, total_time, by = res)

  A_vector <- sapply(
    X = time_vector,
    FUN = island_area,
    total_time = total_time,
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
    num_spec = 1
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
  testit::assert(is.list(out) && is.numeric(unlist(out)))
  testit::assert(length(out) == 2)
  testit::assert(length(out$medians) == 3)
  testit::assert(length(out$means) == 3)
  return(out)
}
