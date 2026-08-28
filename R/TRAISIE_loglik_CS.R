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

    if (is.null(datalist[[1]]$not_present)) {
      loglik <- (datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2) * logp0$loglik
      numimm <- (datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2) + length(datalist) - 1
    } else {
      not_present_by_state <- datalist[[1]]$not_present_by_state
      not_present_NA <- datalist[[1]]$not_present_NA

      # Empirical distribution across OBSERVED states
      p <- not_present_by_state / sum(not_present_by_state)

      # Add the NA species according to the observed-state proportions
      effective_counts_obs <-
        not_present_by_state + not_present_NA * p

      # Expand each observed state over its hidden states
      # assuming hidden states are equally likely
      effective_counts_full <- rep(
        effective_counts_obs / num_hidden_states,
        each = num_hidden_states
      )

      # logp0 has length:
      # num_observed_states * num_hidden_states
      loglik <- sum(effective_counts_full * logp0)

      numimm <-
        sum(not_present_by_state) +
        not_present_NA +
        length(datalist) - 1
    }

    ### pool p0 over the mainland state distribution -> scalar
    w <- effective_counts_full / sum(effective_counts_full)
    logp0_pooled <- log(sum(w * exp(logp0)))

    ### condition on at least one successful colonization
    logcond <- (cond == 1) * log(1 - exp(numimm * logp0_pooled))
    for (i in 2:length(datalist)) {
      datalist[[i]]$type1or2 <- 1
    }
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
