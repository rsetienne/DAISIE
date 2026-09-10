#' testing function, for comparison with DAISIE
#' @description
#' this function calculates the likelihood of observing a singleton endemic species on an island
#' with the trait state `i`, and for which only the estimated maximum and minimum ages of colonization are known.
#' @export
#' @inheritParams default_params_doc
#' @examples
#' library(DAISIE)
#' data("Biwa_datalist")
#' datalist <- Biwa_datalist
#'
#' datalist[[1]]$Mainland_pool_sizes <- c(500, 400)
#' datalist[[1]]$M <- 1000
#'
#' parameter <- list(
#'   c(2.546591, 2.546591, 2.546591, 2.546591),
#'   c(2.678781, 2.678781, 2.678781, 2.678781),
#'   c(0.009326754, 0.009326754, 0.009326754, 0.009326754),
#'   c(1.008583, 1.008583, 1.008583, 1.008583),
#'   matrix(c(
#'     0,    0,    0,  0,
#'     0,    0,    0.00,0.00,
#'     rep(0, 8)
#'   ), nrow = 4, byrow = TRUE),
#'   matrix(1, ncol = 4, nrow = 4)
#' )
#'
#' parameter <- list(
#'   c(2.546591, 1.2, 1, 0.2),
#'   c(2.678781, 2, 1.9, 3),
#'   c(0.009326754, 0.003, 0.002, 0.2),
#'   c(1.008583, 1, 2, 1.5),
#'   matrix(c(
#'     0,    .001,    0.005,  0,
#'     0,    0,    0.002,0.005,
#'     0,    .1000,    0,  0,
#'     0,    0,    0.100,0.00
#'   ), nrow = 4, byrow = TRUE),
#'   1
#' )
#'
#' TRAISIE_logpES_max_min_age_hidden(
#'   datalist                = datalist,
#'   brts                    = c(4, 3.9999, 0.001),
#'   traits                  = 0,
#'   stac                    = 9,
#'   parameter               = parameter,
#'   num_observed_states     = 2,
#'   num_hidden_states       = 2,
#'   atol                    = 1e-15,
#'   rtol                    = 1e-15,
#'   methode                 = "ode45",
#'   trait_mainland_ancestor = c(1, 0),
#'   sampling_fraction       = c(1, 1)
#' )
TRAISIE_logpES_max_min_age_hidden <- function(
    datalist,
    brts,
    parameter,
    traits,
    num_observed_states,
    num_hidden_states,
    trait_mainland_ancestor = NA, # this should contain either a full probability distribution across all states, only the observed states, or NA
    stac,
    sampling_fraction,
    atol = 1e-15,
    rtol = 1e-15,
    methode = "odeint::runge_kutta_cash_karp54"
) {

  calc_Lk_log <- function(i) {
    trait_mainland_ancestor_extended <- rep(0, num_observed_states * num_hidden_states)
    trait_mainland_ancestor_extended[i] <- 1 #set only the trait of interest to 1

    Lk_log <- TRAISIE_logpES_max_min_age_hidden_core(brts,
                                                     parameter               = parameter,
                                                     traits                  = traits,
                                                     num_observed_states     = num_observed_states,
                                                     num_hidden_states       = num_hidden_states,
                                                     trait_mainland_ancestor = trait_mainland_ancestor_extended,
                                                     stac                  = stac,
                                                     sampling_fraction       = sampling_fraction,
                                                     atol                    = atol,
                                                     rtol                    = rtol,
                                                     methode                 = methode)
    return(Lk_log)
  }

  indices_vec <- seq_len(num_observed_states * num_hidden_states)
  Lk_vec <- sapply(indices_vec, calc_Lk_log)

  weights <- TRAISIE_weights(trait_mainland_ancestor,
                             num_observed_states,
                             num_hidden_states,
                             datalist)

  log_Lk <- log(sum(Lk_vec * weights))
  return(list(loglik = log_Lk, lik_states = Lk_vec, weights = weights))
}


#' @keywords internal
TRAISIE_logpES_max_min_age_hidden_core <- function(brts,
                                                   traits,
                                                   stac,
                                                   sampling_fraction = 1,
                                                   parameter,
                                                   trait_mainland_ancestor = NA,
                                                   num_observed_states,
                                                   num_hidden_states,
                                                   atol = 1e-15,
                                                   rtol = 1e-15,
                                                   methode = "odeint::runge_kutta_cash_karp54") {
  t0   <- brts[1]
  tmax <- brts[2]
  tmin <- brts[3]
  tp   <- 0

  # number of unique state
  n <- num_observed_states * num_hidden_states

  #########  interval2 [t_p, tmin]

  m <- length(parameter[[1]])

  initial_conditions2 <- TRAISIE_get_initial_conditions2(stac = stac,
                                                         num_observed_states = num_observed_states,
                                                         num_hidden_states = num_hidden_states,
                                                         trait = traits,
                                                         brts = brts,
                                                         sampling_fraction = sampling_fraction,
                                                         trait_mainland_ancestor = trait_mainland_ancestor)

  # Time sequence for interval [tp, tmin]
  time2 <- c(tp, tmin)

  solution2 <- TRAISIE_solve_branch(interval_func = TRAISIE_interval2,
                                    initial_conditions = initial_conditions2,
                                    time = time2,
                                    parameter = parameter,
                                    methode = methode,
                                    trait_mainland_ancestor = trait_mainland_ancestor,
                                    atol = atol,
                                    rtol = rtol)

  #########  interval3 [tmin, tmax]

  # Initial conditions

  # only use second row, because the first row of solution3 is the initial state
  initial_conditions3_max_min <- c(solution2[2, ][1:n],                                             ### DE: select DE in solution2
                                   rep(0, n),                                                      ### DM1: select DE in solution2
                                   solution2[2, ][(n + 1):(n + n)],                         ### DM2: select DM2 in solution2
                                   solution2[2, ][(n + n + 1):(n + n + n)],                 ### DM3: select DM3 in solution2
                                   solution2[2, ][(n + n + n + 1):(n + n + n + n)],         ### E: select E in solution2
                                   0,                                                              ### DA2
                                   solution2[2, ][length(solution2[2, ])])                           ### DA3: select DA3 in solution2

  initial_conditions3_max_min <- matrix(initial_conditions3_max_min, nrow = 1)

  # Time sequence for interval [tmin, tmax]
  time3 <- c(tmin, tmax)

  solution3 <- TRAISIE_solve_branch(interval_func = TRAISIE_interval3,
                                    initial_conditions = initial_conditions3_max_min,
                                    time = time3,
                                    parameter = parameter,
                                    trait_mainland_ancestor = trait_mainland_ancestor,
                                    methode = methode,
                                    atol = atol,
                                    rtol = rtol)

  #########  interval4 [tmax, t0]

  # Initial conditions

  # only use second row, because the first row of solution3 is the initial state
  initial_conditions4_max_min <- c(solution3[2, ][(n + n + 1):(n + n + n)],                         ### DM1: select DM2 in solution3
                                   solution3[2, ][(n + n + n + n + 1):(n + n + n + n + n)],         ### E: select E in solution3
                                   solution3[2, ][length(solution3[2, ]) - 1])

  initial_conditions4_max_min <- matrix(initial_conditions4_max_min, nrow = 1)

  # Time sequence for interval [tmax, t0]
  time4 <- c(tmax, t0)

  # Solve the system for interval [tmax, t0]
  solution4 <- TRAISIE_solve_branch(interval_func = TRAISIE_interval4,
                                    initial_conditions = initial_conditions4_max_min,
                                    time = time4,
                                    parameter = parameter,
                                    trait_mainland_ancestor = trait_mainland_ancestor,
                                    methode = methode,
                                    atol = atol,
                                    rtol = rtol)

  # Extract log-likelihood
  Lk <- solution4[2, ][length(solution4[2, ])]

  return(Lk)
}
