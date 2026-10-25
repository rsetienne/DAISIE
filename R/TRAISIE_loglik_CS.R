#' testing function, for comparison with DAISIE
#' @description
#' Computes the likelihood of observed island clades and absent mainland lineages.
#' Supports rho-sampling and missing counts per observed trait state.
#' @export
#' @inheritParams default_params_doc
#' @rawNamespace importFrom(RcppParallel, RcppParallelLibs)
TRAISIE_loglik_CS <- function(
    parameter,
    datalist,
    methode = "odeint::runge_kutta_cash_karp54",
    atol = 1e-15,
    rtol = 1e-15,
    num_observed_states,
    num_hidden_states,
    cond = 1,
    num_threads = 1,
    verbose = FALSE,
    sampling = "rho"
) {


  ## will produce error if sampling is not rho or n
  sampling <- match.arg(sampling, c("rho", "n"))




  if (length(datalist) < 2) {
    stop("datalist must contain data and at least one clade.")
  }

  # ------------------------------------------------------------
  # Likelihood for non-present mainland lineages and conditioning
  # ------------------------------------------------------------

  logp0 <- TRAISIE_logp0(
    datalist = datalist,
    parameter = parameter,
    atol = atol,
    rtol = rtol,
    num_observed_states = num_observed_states,
    num_hidden_states = num_hidden_states,
    trait_mainland_ancestor = NA,
    methode = methode
  )
  not_present <- TRAISIE_get_not_present_counts(
    datalist = datalist,
    num_observed_states = num_observed_states,
    num_hidden_states = num_hidden_states
  )
  effective_counts_full <- not_present$effective_counts_full

  # logp0 has length num_observed_states * num_hidden_states
  loglik <- sum(effective_counts_full * logp0)

  numimm <- not_present$total + length(datalist) - 1

  ###
  if (not_present$total > 0) {
    w <- effective_counts_full / sum(effective_counts_full)
  } else {
    # no non-colonising species? then equal weights
    w <- rep(1 / length(logp0), length(logp0))
  }
  logp0_pooled <- log(sum(w * exp(logp0)))

  ### condition on at least one successful colonization
  logcond <- 0
  if (cond %in% c(1, "proper_cond")) {
    logcond <- log(1 - exp(numimm * logp0_pooled))
  }



  loglik <- loglik - logcond

  # ------------------------------------------------------------
  # Likelihood of each observed clade
  # ------------------------------------------------------------

  vec_loglikelihood <- numeric(length(datalist) - 1)

  for (i in seq.int(2, length(datalist))) {

    clade <- datalist[[i]]

    stac <- clade$stac
    brts <- clade$branching_times
    traits <- clade$traits
    phy <- clade$phylogeny
    trait_mainland_ancestor <- clade$root_state

    # Calculate one clade likelihood for the supplied sampling fractions.
    rho_likelihood <- function(sampling_fraction) {
      if (stac %in% c(1, 4)) {
        result <- TRAISIE_logpNE(
          datalist = datalist,
          brts = brts,
          stac = stac,
          traits = traits,
          trait_mainland_ancestor = trait_mainland_ancestor,
          num_observed_states = num_observed_states,
          num_hidden_states = num_hidden_states,
          sampling_fraction = sampling_fraction,
          parameter = parameter,
          atol = atol,
          rtol = rtol,
          methode = methode
        )
      } else if (stac %in% c(2, 3, 5) && length(brts) == 2) {
        result <- TRAISIE_logpES(
          datalist = datalist,
          brts = brts,
          stac = stac,
          traits = traits,
          trait_mainland_ancestor = trait_mainland_ancestor,
          num_observed_states = num_observed_states,
          num_hidden_states = num_hidden_states,
          sampling_fraction = sampling_fraction,
          parameter = parameter,
          atol = atol,
          rtol = rtol,
          methode = methode
        )
      } else if (stac %in% c(2, 3, 5, 6)) {
        result <- TRAISIE_logpEC(
          datalist = datalist,
          brts = brts,
          stac = stac,
          traits = traits,
          trait_mainland_ancestor = trait_mainland_ancestor,
          num_observed_states = num_observed_states,
          num_hidden_states = num_hidden_states,
          sampling_fraction = sampling_fraction,
          parameter = parameter,
          atol = atol,
          rtol = rtol,
          methode = methode,
          phy = phy,
          num_threads = num_threads
        )
      } else if (stac == 8) {
        result <- TRAISIE_logpNE_max_min_age_hidden(
          datalist = datalist,
          brts = brts,
          stac = stac,
          traits = traits,
          trait_mainland_ancestor = trait_mainland_ancestor,
          num_observed_states = num_observed_states,
          num_hidden_states = num_hidden_states,
          sampling_fraction = sampling_fraction,
          parameter = parameter,
          atol = atol,
          rtol = rtol,
          methode = methode
        )
      } else if (stac == 9) {
        result <- TRAISIE_logpES_max_min_age_hidden(
          datalist = datalist,
          brts = brts,
          stac = stac,
          traits = traits,
          trait_mainland_ancestor = trait_mainland_ancestor,
          num_observed_states = num_observed_states,
          num_hidden_states = num_hidden_states,
          sampling_fraction = sampling_fraction,
          parameter = parameter,
          atol = atol,
          rtol = rtol,
          methode = methode
        )
      } else if (stac == 7) {
        stop("stac = 7 needs a likelihood implementation in the underlying TRAISIE code.")
      } else {
        stop("Unknown stac value: ", stac)
      }

      return(result$loglik)
    }

    if (sampling == "rho") {
      vec_loglikelihood[i - 1] <- rho_likelihood(clade$sampling_fraction)
      next
    }



    missing <- clade$missing_species

    if (length(missing) != num_observed_states) {
      stop("Supply one missing-species count per observed trait state.")
    }


    if (all(missing == 0)) {
      vec_loglikelihood[i - 1] <- rho_likelihood(rep(1, num_observed_states))
      next
    }

    if (!stac %in% c(2, 3, 5, 6)) {
      stop("missing-species are only supported for stac 2, 3, 5 and 6")
    }

    # Positive missing counts require complex-valued ODE calculations.
    if (methode != "vode") {
      stop("Positive missing-species counts currently require methode = 'vode'.")
    }

    tip_traits <- as.numeric(traits)


    # Count how many sampled tips are in each observed trait state.
    # Traits are numbered 0, 1, ..., num_observed_states - 1, but tabulate()
    # only counts the numbers 1, 2, 3, ... and silently ignores 0.
    # Adding 1 makes state 0 counted as 1, state 1 counted as 2, etc.,
    # so no tip is lost.
    # Example: traits = c(0, 0, 1, 0)  ->  sampled = c(3, 1)
    #          (3 tips in state 0, 1 tip in state 1)
    sampled <- tabulate(
      tip_traits + 1,
      nbins = num_observed_states
    )

    loglik_rho <- function(rho) {
      rho_likelihood(rho)
    }

    vec_loglikelihood[i - 1] <- TRAISIE_n(
      loglik_rho = loglik_rho,
      sampled = sampled,
      missing = missing,
      radius = pmin(pmax(missing / pmax(sampled + missing, 1), 0.1), 0.95)
    )

  } # End clade loop

  loglik <- loglik + sum(vec_loglikelihood)

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

