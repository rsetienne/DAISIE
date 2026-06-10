#' Testing function for comparison with DAISIE
#'
#' @description
#' This function calculates the likelihood of observing a clade with specified species trait states,
#' given known colonization time. It is designed for comparison with DAISIE-based models.
#'
#' @inheritParams default_params_doc
#'
#' @export
#'
#' @examples
#' library(DAISIE)
#' data("Galapagos_datalist")
#' datalist <- Galapagos_datalist
#' datalist[[1]]$Mainland_pool_sizes <- c(500, 500)
#' datalist[[1]]$M <- 1000
#' i <- 9
#' brts <- datalist[[i]]$branching_times
#' trait <- 0
#' sampling_fraction <- c(1,1)
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
#'   0
#' )
#'
#'
#' parameter <- list(
#'   c(2.546591, 1.2, 1, 0.2),
#'   c(2.678781, 2, 1.9, 3),
#'   c(0.009326754, 0.003, 0.002, 0.2),
#'   c(1.008583, 1, 2, 1.5),
#'   matrix(c(
#'     0,    0.1,    0.05,  0,
#'     0.33,    0,    0.000,0.0086,
#'     0.005,    000,    0,  0.005,
#'     0,   0.5,  0.35,0.00
#'   ), nrow = 4, byrow = TRUE),
#'   1
#' )
#'
#' status = 2
#' TRAISIE_logpES(
#'   datalist              = datalist,
#'   brts                    = brts,
#'   trait                   = trait,
#'   status                  = status,
#'   sampling_fraction       = sampling_fraction,
#'   parameter               = parameter,
#'   trait_mainland_ancestor = c(1,0),
#'   num_observed_states     = 2,
#'   num_hidden_states       = 2,
#'   atol                    = 1e-15,
#'   rtol                    = 1e-15,
#'   methode                 = "odeint::runge_kutta_cash_karp54")



TRAISIE_logpES <- function(
    datalist,
    brts,
    parameter,
    trait,
    num_observed_states,
    num_hidden_states,
    trait_mainland_ancestor = NA, #this should contain either a full probability distribution across all states, only the observed states, or NA
    status,
    sampling_fraction,
    atol = 1e-15,
    rtol = 1e-15,
    methode = "odeint::runge_kutta_cash_karp54"
) {

  calc_Lk_log <- function(i) {
    trait_mainland_ancestor_extended <- rep(0, num_observed_states * num_hidden_states)
    trait_mainland_ancestor_extended[i] <- 1 #set only the trait of interest to 1

    Lk_log <-  TRAISIE_logpES_core(brts                    = brts,
                                   parameter               = parameter,
                                   trait                   = trait,
                                   num_observed_states     = num_observed_states,
                                   num_hidden_states       = num_hidden_states,
                                   trait_mainland_ancestor = trait_mainland_ancestor_extended,
                                   status                  = status,
                                   sampling_fraction       = sampling_fraction,
                                   atol                    = atol,
                                   rtol                    = rtol,
                                   methode                 = methode,
                                   use_Rcpp                = use_Rcpp)
    return(Lk_log)
  }

  indices_vec <- seq_len(num_observed_states * num_hidden_states)
  Lk_vec <- sapply(indices_vec, calc_Lk_log)

  ## added !all(is.na(trait_mainland_ancestor)) because when trait_mainland_ancestor = NA,  length(trait_mainland_ancestor) = length(trait_mainland_ancestor_extended) = 1
  if (!all(is.na(trait_mainland_ancestor)) && length(trait_mainland_ancestor) == num_observed_states * num_hidden_states) { #this is the case where a full probability distribution is specified across all observed and hidden states

    weights <- trait_mainland_ancestor / sum(trait_mainland_ancestor)
  }  else {

    if (all(is.numeric(trait_mainland_ancestor))) { # this is the case when only a probability distribution is specified for the observed states; this could be c(M0/M, M1/M)

      s <- numeric(num_observed_states * num_hidden_states)
      # you could also do s <- c() and use line 92

      weights <- c()
      for (j in 1:length(trait_mainland_ancestor)) {
        s[((j - 1) * num_hidden_states + 1):(j * num_hidden_states)] <- rep(trait_mainland_ancestor[j], num_hidden_states)


      }
      weights <- s / sum(s)

    }else { # this is the case where nothing is provided, i.e. NA
      Mp <- datalist[[1]]$Mainland_pool_sizes
      M <-  datalist[[1]]$M
      num_hidden_states <- num_hidden_states
      weights <- TRAISIE_compute_mainland_weights(Mp, M, num_hidden_states)

    }
  }
  log_Lk <- log(sum(Lk_vec * weights))
  return(list(loglik = log_Lk, lik_states = Lk_vec, weights = weights))
}


#' @keywords internal
TRAISIE_logpES_core <- function(brts,
                                status,
                                trait,
                                sampling_fraction,
                                num_observed_states,
                                num_hidden_states,
                                trait_mainland_ancestor = NA,
                                parameter,
                                atol  = 1e-15,
                                rtol  = 1e-15,
                                methode = "odeint::runge_kutta_cash_karp54",
                                use_Rcpp = 0) {


  TRAISIE_check_arguments(brts, parameter,
                          phy = 0,
                          trait,
                          num_observed_states,
                          num_hidden_states,
                          status,
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

  # Number of states in the system

  # Solve for interval [tp, t2] (stem phase)


  # Run appropriate sequence of intervals
  if ((status == 2 || status == 3) && length(brts) == 2) {
    initial_conditions2 <- TRAISIE_get_initial_conditions2(status = status,
                                                           num_observed_states = num_observed_states,
                                                           num_hidden_states = num_hidden_states,
                                                           trait = trait,
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
                                      rtol = rtol,
                                      use_Rcpp = use_Rcpp)


    initial_conditions4 <- TRAISIE_get_initial_conditions4(status = status,
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
                                      rtol = rtol,
                                      use_Rcpp = use_Rcpp)
  }

  if (status == 5) {
    initial_conditions3 <- TRAISIE_get_initial_conditions3(status = status,
                                                           num_observed_states = num_observed_states,
                                                           num_hidden_states = num_hidden_states,
                                                           trait = trait,
                                                           sampling_fraction = sampling_fraction)
    solution3 <- TRAISIE_solve_branch(interval_func = TRAISIE_interval3,
                                      initial_conditions = initial_conditions3,
                                      time = time3,
                                      parameter = parameter,
                                      trait_mainland_ancestor = trait_mainland_ancestor,
                                      methode = methode,
                                      atol = atol,
                                      rtol = rtol,
                                      use_Rcpp = use_Rcpp)

    initial_conditions4 <- TRAISIE_get_initial_conditions4(status = status,
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
                                      rtol = rtol,
                                      use_Rcpp = use_Rcpp)
  }

  # Extract log-likelihood from final solution
  Lk <- solution4[2, length(solution4[2, ])]

  return(Lk)
}
