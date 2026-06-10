#' testing function, for comparison with DAISIE
#' @description
#' This function compute the likelihood that all species that colonize the island
#' have gone extinct prior to the present.
#' @export
#' @inheritParams default_params_doc
#' @examples
#' #Load DAISIE package and data
#' library(DAISIE)
#' data("Galapagos_datalist")
#' datalist <- Galapagos_datalist
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
#'   ), nrow = 4),
#'   1
#' )
#' datalist[[1]]$Mainland_pool_sizes <- c(500, 400)
#' datalist[[1]]$M <- 1000
#'
#'
#' TRAISIE_logp0(
#'   datalist,
#'   parameter               = parameter,
#'   trait_mainland_ancestor = NA,
#'   num_observed_states     = 2,
#'   num_hidden_states       = 2,
#'   atol                    = 1e-15,
#'   rtol                    = 1e-15,
#'   methode                 = "odeint::runge_kutta_cash_karp54",
#'   use_Rcpp                = 2)
#'
TRAISIE_logp0 <- function(
    datalist,
    parameter,
    atol = 1e-15,
    rtol = 1e-15,
    num_observed_states,
    num_hidden_states,
    trait_mainland_ancestor = NA,
    methode = "odeint::runge_kutta_cash_karp54",
    use_Rcpp = 2) {

  calc_Lk_log <- function(i) {
    trait_mainland_ancestor_extended <- rep(0, num_observed_states * num_hidden_states)
    trait_mainland_ancestor_extended[i] <- 1 #set only the trait of interest to 1

    Lk_log <- TRAISIE_logp0_core(datalist,
                                 parameter,
                                 atol = 1e-15,
                                 rtol = 1e-15,
                                 num_observed_states,
                                 num_hidden_states,
                                 trait_mainland_ancestor = trait_mainland_ancestor_extended,
                                 methode = methode,
                                 use_Rcpp = use_Rcpp)
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
  return(log_Lk)
}


#' @keywords internal
TRAISIE_logp0_core <- function(datalist,
                               parameter,
                               atol = 1e-15,
                               rtol = 1e-15,
                               num_observed_states,
                               num_hidden_states,
                               trait_mainland_ancestor = NA,
                               methode = "odeint::runge_kutta_cash_karp54",
                               use_Rcpp = 0) {

  n <- num_observed_states * num_hidden_states
  t0 <- datalist[[1]]$island_age
  tp <- 0

  ######### interval4 [t_p, t_0]

  initial_conditions40 <- c(rep(0, n),  ### DM1
                            rep(0, n),  ### E
                            1)          ### DA1

  # Time sequence for interval [tp, t0]
  time4 <- c(tp, t0)

  # Solve the system for interval [tp, t1]
  solution4 <- TRAISIE_solve_branch(interval_func = TRAISIE_interval4,
                                    initial_conditions = initial_conditions40,
                                    time = time4,
                                    parameter = parameter,
                                    trait_mainland_ancestor = trait_mainland_ancestor,
                                    methode = methode,
                                    atol = atol,
                                    rtol = rtol,
                                    use_Rcpp = use_Rcpp)

  # Extract log-likelihood
  Lk <- solution4[2, ][length(solution4[2, ])]

  return(Lk)
}



TRAISIE_logp0_core <- function(datalist,
                               parameter,
                               atol = 1e-15,
                               rtol = 1e-15,
                               num_observed_states,
                               num_hidden_states,
                               trait_mainland_ancestor = NA,
                               methode = "odeint::runge_kutta_cash_karp54",
                               use_Rcpp = 0) {

  n <- num_observed_states * num_hidden_states
  t0 <- datalist[[1]]$island_age
  tp <- 0

  #########  interval4 [t_p, t_0]

  initial_conditions40 <- c(rep(0, n),  ### DM1
                            rep(0, n),  ### E
                            1)          ### DA1

  # Time sequence for interval [tp, t0]
  time4 <- c(tp, t0)

  # Solve the system for interval [tp, t1]
  solution4 <- TRAISIE_solve_branch(interval_func = TRAISIE_interval4,
                                    initial_conditions = initial_conditions40,
                                    time = time4,
                                    parameter = parameter,
                                    trait_mainland_ancestor = trait_mainland_ancestor,
                                    methode = methode,
                                    atol = atol,
                                    rtol = rtol,
                                    use_Rcpp = use_Rcpp)

  # Extract log-likelihood
  Lk <- solution4[2, ][length(solution4[2, ])]

  return(Lk)
}
