#' testing function, for comparison with DAISIE
#' @description
#' This function compute the likelihood that all species that colonize the island
#' have gone extinct prior to the present.
#' @export
#' @inheritParams default_params_doc
#' @rawNamespace importFrom(RcppParallel, RcppParallelLibs)
TRAISIE_loglik_CS <- function(parameter,
                              datalist,
                              methode = "odeint::runge_kutta_cash_karp54",
                              atol = 1e-15,
                              rtol = 1e-15,
                              num_observed_states,
                              num_hidden_states,
                              cond = 1,
                              num_threads = 1,
                              verbose = FALSE) {
  logcond <- 0 # default value gives no effect

  if (length(parameter) >= 6) {
    logp0 <- TRAISIE_logp0(datalist = datalist,
                           parameter = parameter,
                           atol = atol,
                           rtol = rtol,
                           num_observed_states = num_observed_states,
                           num_hidden_states = num_hidden_states,
                           trait_mainland_ancestor =  NA,
                           methode = methode)

    not_present <- TRAISIE_get_not_present_counts(
      datalist = datalist,
      num_observed_states = num_observed_states,
      num_hidden_states = num_hidden_states
    )
    effective_counts_full <- not_present$effective_counts_full

    # logp0 has length num_observed_states * num_hidden_states
    loglik <- sum(effective_counts_full * logp0)

    numimm <- not_present$total + length(datalist) - 1

    ### pool p0 over the mainland state distribution -> scalar
    if (not_present$total > 0) {
      w <- effective_counts_full / sum(effective_counts_full)
    } else {
      # no non-colonising species: fall back to equal weights
      w <- rep(1 / length(logp0), length(logp0))
    }
    logp0_pooled <- log(sum(w * exp(logp0)))

    ### condition on at least one successful colonization
    logcond <- 0
    if (cond %in% c(1, "proper_cond")) {
      logcond <- log(1 - exp(numimm * logp0_pooled))
    }

  }

  if (!is.list(parameter) || length(parameter) < 6) {
    stop("'parameter' must be a list with 6 elements ",
         "(lambda_c, mu, gamma, lambda_a, q, p); got ", ...)
  }

  loglik <- loglik - logcond

  vec_loglikelihood <- rep(NA, length(datalist) - 1) # first entry is not data
  for (i in 2:length(datalist)) {
    stac <- datalist[[i]]$stac
    brts <- datalist[[i]]$branching_times
    traits <- datalist[[i]]$traits

    sampling_fraction <- datalist[[i]]$sampling_fraction

    phy <- datalist[[i]]$phylogeny

    trait_mainland_ancestor <- datalist[[i]]$root_state

    if (stac %in% c(1, 4)) {
      loglikelihood <- TRAISIE_logpNE(datalist = datalist,
                                      brts = brts,
                                      stac = stac,
                                      traits = traits,
                                      trait_mainland_ancestor = trait_mainland_ancestor,
                                      num_observed_states = num_observed_states,
                                      num_hidden_states = num_hidden_states,
                                      sampling_fraction = sampling_fraction,
                                      parameter = parameter,
                                      atol  = atol,
                                      rtol  = rtol,
                                      methode = methode)
    } else if (stac %in% c(2, 5)) {
      if (length(brts) == 2) {
        loglikelihood <- TRAISIE_logpES(datalist = datalist,
                                        brts = brts,
                                        stac = stac,
                                        traits = traits,
                                        sampling_fraction = sampling_fraction,
                                        trait_mainland_ancestor = trait_mainland_ancestor,
                                        num_observed_states = num_observed_states,
                                        num_hidden_states = num_hidden_states,
                                        parameter = parameter,
                                        atol  = atol,
                                        rtol  = rtol,
                                        methode = methode)
      } else {
        loglikelihood <- TRAISIE_logpEC(datalist = datalist,
                                        brts = brts,
                                        parameter = parameter,
                                        phy = phy,
                                        traits = traits,
                                        num_observed_states = num_observed_states,
                                        num_hidden_states = num_hidden_states,
                                        trait_mainland_ancestor = trait_mainland_ancestor,
                                        stac = stac,
                                        sampling_fraction = sampling_fraction,
                                        atol  = atol,
                                        rtol  = rtol,
                                        methode = methode,
                                        num_threads = num_threads)
      }
    } else if (stac == 3) {
      if (length(brts) == 2) {
        loglikelihood <- TRAISIE_logpES(datalist = datalist,
                                        brts = brts,
                                        stac = stac,
                                        traits = traits,
                                        sampling_fraction = sampling_fraction,
                                        trait_mainland_ancestor = trait_mainland_ancestor,
                                        num_observed_states = num_observed_states,
                                        num_hidden_states = num_hidden_states,
                                        parameter = parameter,
                                        atol  = atol,
                                        rtol  = rtol,
                                        methode = methode)
      } else {
        loglikelihood <- TRAISIE_logpEC(datalist = datalist,
                                        brts = brts,
                                        parameter = parameter,
                                        phy = phy,
                                        traits = traits,
                                        num_observed_states = num_observed_states,
                                        num_hidden_states = num_hidden_states,
                                        trait_mainland_ancestor = trait_mainland_ancestor,
                                        stac = stac,
                                        sampling_fraction = sampling_fraction,
                                        atol  = atol,
                                        rtol  = rtol,
                                        methode = methode)
      }
    } else if (stac == 6) {
      loglikelihood <- TRAISIE_logpEC(datalist = datalist,
                                      brts = brts,
                                      parameter = parameter,
                                      phy = phy,
                                      traits = traits,
                                      num_observed_states = num_observed_states,
                                      num_hidden_states = num_hidden_states,
                                      trait_mainland_ancestor = trait_mainland_ancestor,
                                      stac = stac,
                                      sampling_fraction = sampling_fraction,
                                      atol  = atol,
                                      rtol  = rtol,
                                      methode = methode,
                                      num_threads = num_threads)
    } else if (stac == 8) {
      loglikelihood <-
        TRAISIE_logpNE_max_min_age_hidden(datalist = datalist,
                                          brts = brts,
                                          traits = traits,
                                          stac = stac,
                                          parameter = parameter,
                                          trait_mainland_ancestor = trait_mainland_ancestor,
                                          num_observed_states = num_observed_states,
                                          num_hidden_states = num_hidden_states,
                                          sampling_fraction = sampling_fraction,
                                          atol  = atol,
                                          rtol  = rtol,
                                          methode = methode)
    } else if (stac == 9) {
      loglikelihood <-
        TRAISIE_logpES_max_min_age_hidden(datalist = datalist,
                                          brts = brts,
                                          traits = traits,
                                          sampling_fraction = sampling_fraction,
                                          stac = stac,
                                          parameter = parameter,
                                          trait_mainland_ancestor = trait_mainland_ancestor,
                                          num_observed_states = num_observed_states,
                                          num_hidden_states = num_hidden_states,
                                          atol  = atol,
                                          rtol  = rtol,
                                          methode = methode)
    } else {
      stop("Unknown stac value: ", stac)
    }

    vec_loglikelihood[i - 1] <- loglikelihood$loglik
  }

  loglik <- sum(vec_loglikelihood) + loglik
  return(loglik)

}



#' Effective number of non-colonising mainland species per (observed x hidden)
#' state
#'
#' Reads `not_present_by_state` and `not_present_NA` from `datalist[[1]]`,
#' checks them, distributes the species with unknown state over the observed
#' states in proportion to the known counts, and splits each observed state
#' equally over its hidden states.
#' @return list with `effective_counts_full` (length
#'   num_observed_states * num_hidden_states) and `total`.
#' @keywords internal
#' @noRd
TRAISIE_get_not_present_counts <- function(datalist,
                                           num_observed_states,
                                           num_hidden_states) {
  by_state <- datalist[[1]]$not_present_by_state
  n_na     <- datalist[[1]]$not_present_NA

  if (is.null(by_state)) {
    msg <- "datalist[[1]]$not_present_by_state is missing."
    if (!is.null(datalist[[1]]$not_present) ||
        !is.null(datalist[[1]]$not_present_type1)) {
      msg <- paste(msg, "This looks like a DAISIE datalist (not_present /",
                   "not_present_type1/2); TRAISIE needs the number of",
                   "non-colonising mainland species per trait state.")
    }
    stop(msg, call. = FALSE)
  }
  if (is.null(n_na)) n_na <- 0   # no species with unknown trait state

  if (!is.numeric(by_state) || length(by_state) != num_observed_states) {
    stop("datalist[[1]]$not_present_by_state must be numeric with one entry ",
         "per observed state (", num_observed_states, "); got length ",
         length(by_state), ".", call. = FALSE)
  }
  if (!is.numeric(n_na) || length(n_na) != 1) {
    stop("datalist[[1]]$not_present_NA must be a single number.",
         call. = FALSE)
  }
  if (anyNA(by_state) || is.na(n_na) || any(by_state < 0) || n_na < 0) {
    stop("not_present_by_state and not_present_NA must be non-negative ",
         "and not NA.", call. = FALSE)
  }

  # distribution over OBSERVED states; uniform if no species have a known state
  p <- if (sum(by_state) > 0) by_state / sum(by_state) else
    rep(1 / num_observed_states, num_observed_states)

  effective_counts_obs  <- by_state + n_na * p
  effective_counts_full <- rep(effective_counts_obs / num_hidden_states,
                               each = num_hidden_states)

  list(effective_counts_full = effective_counts_full,
       total = sum(by_state) + n_na)
}
