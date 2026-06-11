
#' This function calculates the likelihood of observing a non-endemic lineage with specified species trait states,
#'
#'
#' @inheritParams default_params_doc
#'
#' @export
#'
#' @examples
#' library(DAISIE)
#' data("Galapagos_datalist")
#' datalist <- Galapagos_datalist
#' datalist[[1]]$M0 <- 500
#' datalist[[1]]$M1 <- 400
#' i <- 3
#' brts <- datalist[[i]]$branching_times
#' traits <- 0
#'
#' parameter <- list(
#'   c(2.546591, 1.2, 1, 0.2),
#'   c(2.678781, 2, 1.9, 3),
#'   c(0.009326754, 0.003, 0.002, 0.2),
#'   c(1.008583, 1, 2, 1.5),
#'   matrix(c(
#'     0,    .001,    0.005,  0,
#'     .001,    0,    0.000,0.005,
#'     0.005,    000,    0,  0.005,
#'     0,   0.005,  0.005,0.00
#'   ), nrow = 4, byrow = TRUE),
#'   1
#' )
#' brts <- datalist[[i]]$branching_times
#' TRAISIE_logpNE(
#'   datalist                = datalist,
#'   brts                    = brts,
#'   traits                  = traits,
#'   stac                    = 4,
#'   parameter               = parameter,
#'   sampling_fraction       = c(1,1),
#'   trait_mainland_ancestor = c(1,0),
#'   num_observed_states     = 2,
#'   num_hidden_states       = 2,
#'   atol                    = 1e-15,
#'   rtol                    = 1e-15,
#'   methode                 = "odeint::runge_kutta_cash_karp54")
#'
#'
TRAISIE_logpNE <- function(
    datalist,
    brts,
    parameter,
    traits,
    num_observed_states,
    num_hidden_states,
    trait_mainland_ancestor = NA, #this should contain either a full probability distribution across all states, only the observed states, or NA
    stac,
    sampling_fraction,
    atol = 1e-15,
    rtol = 1e-15,
    methode = "odeint::runge_kutta_cash_karp54"
) {
  calc_Lk_log <- function(i) {
    trait_mainland_ancestor_extended <- rep(0, num_observed_states * num_hidden_states)
    trait_mainland_ancestor_extended[i] <- 1 #set only the trait of interest to 1

    Lk_log <-  TRAISIE_logpNE_core(brts                    = brts,
                                   parameter               = parameter,
                                   traits                  = traits,
                                   num_observed_states     = num_observed_states,
                                   num_hidden_states       = num_hidden_states,
                                   trait_mainland_ancestor = trait_mainland_ancestor_extended,
                                   stac                    = stac,
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
TRAISIE_logpNE_core <- function(brts,
                                stac,
                                traits,
                                sampling_fraction,
                                num_observed_states,
                                num_hidden_states,
                                trait_mainland_ancestor = NA,
                                parameter,
                                atol  = 1e-15,
                                rtol  = 1e-15,
                                methode = "odeint::runge_kutta_cash_karp54") {

  TRAISIE_check_arguments(brts = brts,
                          parameter = parameter,
                          phy = 1,
                          traits = traits,
                          num_observed_states = num_observed_states,
                          num_hidden_states = num_hidden_states,
                          stac = stac,
                          sampling_fraction = sampling_fraction)
  # Unpack times from brts
  t0   <- brts[1]
  tmax <- brts[2]
  t1   <- brts[2]

  tp   <- 0

  # Time intervals

  time2 <- c(tp, t1)
  time3 <- c(tp, tmax)
  time4 <- c(tmax, t0)

  # Solve for interval [tp, t2] (stem phase)

  # Run appropriate sequence of intervals
  if (stac == 4) {
    initial_conditions2 <- TRAISIE_get_initial_conditions2(stac = stac,
                                                           trait = traits,
                                                           num_observed_states = num_observed_states,
                                                           num_hidden_states = num_hidden_states,
                                                           brts = brts,
                                                           sampling_fraction = sampling_fraction,
                                                           trait_mainland_ancestor = trait_mainland_ancestor)
    solution2 <- TRAISIE_solve_branch(interval_func = TRAISIE_interval2,
                                      initial_conditions = initial_conditions2,
                                      time = time2,
                                      parameter = parameter,
                                      trait_mainland_ancestor = trait_mainland_ancestor,
                                      methode = methode,
                                      atol = atol,
                                      rtol = rtol)

    initial_conditions4 <- TRAISIE_get_initial_conditions4(stac = stac,
                                                           solution = solution2,
                                                           parameter = parameter,
                                                           trait_mainland_ancestor = trait_mainland_ancestor,
                                                           num_observed_states = num_observed_states,
                                                           num_hidden_states = num_hidden_states)
    solution4 <- TRAISIE_solve_branch(interval_func = TRAISIE_interval4,
                                      initial_conditions = initial_conditions4,
                                      time = time4,
                                      parameter = parameter,
                                      trait_mainland_ancestor = trait_mainland_ancestor,
                                      methode = methode,
                                      atol = atol,
                                      rtol = rtol)
  }

  if (stac == 1) {
    initial_conditions3 <- TRAISIE_get_initial_conditions3(stac = stac,
                                                           num_observed_states = num_observed_states,
                                                           num_hidden_states = num_hidden_states,
                                                           trait = traits,
                                                           sampling_fraction = sampling_fraction)
    solution3 <- TRAISIE_solve_branch(interval_func = TRAISIE_interval3,
                                      initial_conditions = initial_conditions3,
                                      time = time3,
                                      parameter = parameter,
                                      trait_mainland_ancestor = trait_mainland_ancestor,
                                      methode = methode,
                                      atol = atol,
                                      rtol = rtol)

    initial_conditions4 <- TRAISIE_get_initial_conditions4(stac = stac,
                                                           solution = solution3,
                                                           parameter = parameter,
                                                           trait_mainland_ancestor = trait_mainland_ancestor,
                                                           num_observed_states = num_observed_states,
                                                           num_hidden_states = num_hidden_states)
    solution4 <- TRAISIE_solve_branch(interval_func = TRAISIE_interval4,
                                      initial_conditions = initial_conditions4,
                                      time = time4,
                                      parameter = parameter,
                                      trait_mainland_ancestor = trait_mainland_ancestor,
                                      methode = methode,
                                      atol = atol,
                                      rtol = rtol)
  }

  # Extract log-likelihood from final solution
  Lk <- solution4[2, length(solution4[2, ])]
  return(Lk)
}
