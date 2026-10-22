#' testing function, for comparison with DAISIE
#' @description
#' This function compute the likelihood that all species that colonize the island
#' have gone extinct prior to the present.
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

  sampling <- match.arg(sampling, c("rho", "n"))

  if (sampling == "n" && methode != "vode") {
    stop("n-sampling currently requires methode = 'vode'.")
  }

  if (length(parameter) < 6L) {
    stop("parameter must contain the six TRAISIE parameter components.")
  }

  if (length(datalist) < 2L) {
    stop("datalist must contain metadata and at least one clade.")
  }

  # ------------------------------------------------------------
  # Absent mainland lineages and conditioning
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

  not_present_by_state <- datalist[[1]]$not_present_by_state
  not_present_NA <- datalist[[1]]$not_present_NA

  # Preserve the existing allocation of unknown mainland states.
  p <- not_present_by_state / sum(not_present_by_state)

  effective_counts_obs <-
    not_present_by_state + not_present_NA * p

  effective_counts_full <- rep(
    effective_counts_obs / num_hidden_states,
    each = num_hidden_states
  )

  loglik <- sum(effective_counts_full * logp0)

  numimm <-
    sum(not_present_by_state) +
    not_present_NA +
    length(datalist) - 1L

  logcond <- 0

  if (cond %in% c(1, "proper_cond")) {
    w <- effective_counts_full / sum(effective_counts_full)
    logp0_pooled <- log(sum(w * exp(logp0)))

    logcond <- log(1 - exp(numimm * logp0_pooled))
  }

  loglik <- loglik - logcond

  # ------------------------------------------------------------
  # Likelihood of each observed clade
  # ------------------------------------------------------------

  vec_loglikelihood <- numeric(length(datalist) - 1L)

  for (i in seq.int(2L, length(datalist))) {

    clade <- datalist[[i]]

    stac <- clade$stac
    brts <- clade$branching_times
    traits <- clade$traits
    phy <- clade$phylogeny
    trait_mainland_ancestor <- clade$root_state

    # Evaluate the existing rho-sampling clade likelihood.
    rho_likelihood <- function(sampling_fraction) {

      args <- list(
        datalist = datalist,
        brts = brts,
        stac = stac,
        traits = traits,
        trait_mainland_ancestor = trait_mainland_ancestor,
        num_observed_states = num_observed_states,
        num_hidden_states = num_hidden_states,
        sampling_fraction = as.vector(sampling_fraction),
        parameter = parameter,
        atol = atol,
        rtol = rtol,
        methode = methode
      )

      if (stac %in% c(1, 4)) {

        likelihood_function <- TRAISIE_logpNE

      } else if (stac %in% c(2, 3, 5)) {

        if (length(brts) == 2L) {
          likelihood_function <- TRAISIE_logpES
        } else {
          likelihood_function <- TRAISIE_logpEC
          args$phy <- phy
          args$num_threads <- num_threads
        }

      } else if (stac == 6) {

        likelihood_function <- TRAISIE_logpEC
        args$phy <- phy
        args$num_threads <- num_threads

      } else if (stac == 8) {

        likelihood_function <- TRAISIE_logpNE_max_min_age_hidden

      } else if (stac == 9) {

        likelihood_function <- TRAISIE_logpES_max_min_age_hidden

      } else {

        stop("Unknown stac value: ", stac)
      }

      result <- do.call(likelihood_function, args)

      result$loglik
    }

    # ----------------------------------------------------------
    # Rho-sampling
    # ----------------------------------------------------------

    if (sampling == "rho") {

      vec_loglikelihood[i - 1L] <- rho_likelihood(
        clade$sampling_fraction
      )

      next
    }

    # ----------------------------------------------------------
    # n-sampling
    # ----------------------------------------------------------

    missing <- clade$missing_species

    if (is.null(missing)) {
      stop("Clade ", i, " requires missing_species.")
    }

    missing <- as.numeric(missing)

    if (any(!is.finite(missing)) ||
        any(missing < 0 | missing != floor(missing))) {
      stop(
        "Clade ", i,
        ": missing_species must contain nonnegative integers."
      )
    }

    if (!length(missing) %in% c(1L, num_observed_states)) {
      stop(
        "Clade ", i,
        ": supply a total missing count or one count ",
        "per observed trait state."
      )
    }

    # Complete sampling needs no coefficient extraction.
    if (all(missing == 0)) {

      vec_loglikelihood[i - 1L] <- rho_likelihood(
        rep(1, num_observed_states)
      )

      next
    }

    # Current scope for positive missing counts.
    if (!stac %in% c(2, 5, 6)) {
      stop(
        "Clade ", i,
        ": positive missing counts are not yet supported ",
        "for stac = ", stac
      )
    }

    tip_traits <- as.numeric(traits)

    number_sampled <- if (inherits(phy, "phylo")) {
      length(phy$tip.label)
    } else {
      length(tip_traits)
    }

    if (number_sampled < 1L ||
        length(tip_traits) != number_sampled) {
      stop(
        "Clade ", i,
        ": provide one trait entry per sampled species."
      )
    }

    known_traits <- tip_traits[!is.na(tip_traits)]

    if (any(!is.finite(known_traits)) ||
        any(known_traits != floor(known_traits)) ||
        any(known_traits < 0 |
            known_traits >= num_observed_states)) {
      stop("Clade ", i, ": invalid observed trait states.")
    }

    if (length(missing) == 1L) {

      # Total number missing; their observed states are unspecified.
      sampled <- number_sampled

      loglik_rho <- function(rho) {
        rho_likelihood(rep(rho, num_observed_states))
      }

    } else {

      # Missing counts supplied separately by observed trait state.
      if (anyNA(tip_traits)) {
        stop(
          "Clade ", i,
          ": trait-specific n-sampling currently requires ",
          "known tip traits. Use a scalar missing_species ",
          "count for total-count sampling."
        )
      }

      sampled <- tabulate(
        as.integer(tip_traits) + 1L,
        nbins = num_observed_states
      )

      loglik_rho <- function(rho) {
        rho_likelihood(rho)
      }
    }

    vec_loglikelihood[i - 1L] <- TRAISIE_n(
      loglik_rho = loglik_rho,
      sampled = sampled,
      missing = missing,
      radius = 0.1
    )

  } # End clade loop

  loglik <- loglik + sum(vec_loglikelihood)

  return(loglik)
}
