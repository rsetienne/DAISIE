test_that("DAISIE_DE and DAISIE agree", {
  skip_on_cran()
  data(made_up_datalist)
  datalist <- made_up_datalist

  # =============================================================================
  # Create all new entries and add them to datalist
  # =============================================================================

  make_artificial_entry <- function(
    missing_species,
    stac = 3,
    branching_times = c(15.80, 4.6934329)
  ) {
    list(
      branching_times = branching_times,
      stac = stac,
      missing_species = missing_species
    )
  }

  modify_entry <- function(entry, missing_species = NULL, stac = NULL) {

    new_entry <- entry

    if (!is.null(missing_species)) {
      new_entry$missing_species <- missing_species
    }

    if (!is.null(stac)) {
      new_entry$stac <- stac
    }

    return(new_entry)
  }

  # =============================================================================
  # Create all entries
  # =============================================================================

  new_entries <- list(

    # ---------------------------------------------------------------------------
    # Artificial stac 3 entries with 2 branching times
    # ---------------------------------------------------------------------------

    make_artificial_entry(
      missing_species = 0,
      stac = 3
    ),

    make_artificial_entry(
      missing_species = 2,
      stac = 3
    ),

    # ---------------------------------------------------------------------------
    # Artificial stac 7 entries with 2 branching times
    # ---------------------------------------------------------------------------

    make_artificial_entry(
      missing_species = 0,
      stac = 7
    ),

    make_artificial_entry(
      missing_species = 2,
      stac = 7
    ),

    # ---------------------------------------------------------------------------
    # Existing datalist-based entries
    # ---------------------------------------------------------------------------

    # Same as datalist[[4]], but missing_species = 2
    modify_entry(
      entry = datalist[[4]],
      missing_species = 2
    ),

    # Same as datalist[[9]], but missing_species = 2 and stac = 5
    modify_entry(
      entry = datalist[[9]],
      missing_species = 2,
      stac = 5
    ),

    # Same as datalist[[6]], but missing_species = 2 and stac = 6
    modify_entry(
      entry = datalist[[6]],
      missing_species = 2,
      stac = 6
    ),

    # Same as datalist[[6]], but stac = 6
    modify_entry(
      entry = datalist[[6]],
      stac = 6
    ),

    # ---------------------------------------------------------------------------
    # Artificial stac 3 entries with more than 2 branching times
    # ---------------------------------------------------------------------------

    make_artificial_entry(
      missing_species = 0,
      stac = 3,
      branching_times = datalist[[6]]$branching_times
    ),

    make_artificial_entry(
      missing_species = 2,
      stac = 3,
      branching_times = datalist[[6]]$branching_times
    ),

    # ---------------------------------------------------------------------------
    # Artificial stac 7 entries with more than 2 branching times
    # ---------------------------------------------------------------------------

    make_artificial_entry(
      missing_species = 0,
      stac = 7,
      branching_times = c(15.8, 2, 1)
    ),

    make_artificial_entry(
      missing_species = 2,
      stac = 7,
      branching_times = c(15.8, 2, 1)
    )
  )

  # =============================================================================
  # Add all entries to datalist at once
  # =============================================================================

  datalist <- c(datalist, new_entries)

  pars1 <- c(2.546591, 2.678781, 2.678781, 0.009326754, 1.008583)

  ll1 <- DAISIE:::DAISIE_loglik_CS(
    pars1 = pars1,
    pars2 = c(100, 11, 0, 2),
    datalist = datalist,
    methode = "odeint::runge_kutta_cash_karp54",
    abstolint = 1e-16,
    reltolint = 1e-16,
    CS_version = list(
      model = 1,
      function_to_optimize = "DAISIE_DE",
      sampling = "n"
    )
  )

  pars1[3] <- Inf

  ll2 <- DAISIE:::DAISIE_loglik_CS(
    pars1 = pars1,
    pars2 = c(100, 1, 0, 2),
    datalist = datalist,
    methode = "odeint::runge_kutta_cash_karp54",
    abstolint = 1e-16,
    reltolint = 1e-16,
    CS_version = list(
      model = 1,
      function_to_optimize = "DAISIE",
      sampling = "n"
    )
  )

  testthat::expect_equal(ll1, ll2, tolerance = 1e-4)
})
