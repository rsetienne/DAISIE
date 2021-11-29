test_that("silent with empty island with correct output", {
  pars <- c(0.4, 0.2, 10, 0.0001, 0.5)
  time <- 1
  mainland_n <- 1
  verbose <- FALSE
  sample_freq <- 1
  set.seed(1)
  island_replicates <- list()
  out <- list()
  out[[1]] <- DAISIE:::DAISIE_sim_core_cr(
    time = time,
    pars = pars,
    mainland_n = mainland_n,
    nonoceanic_pars = c(0, 0),
    hyper_pars = create_hyper_pars(d = 0, x = 0),
    area_pars = create_area_pars(
      max_area = 1,
      current_area = 1,
      proportional_peak_t = 0,
      total_island_age = 0,
      sea_level_amplitude = 0,
      sea_level_frequency = 0,
      island_gradient_angle = 0
    )
  )
  island_replicates[[1]] <- out
  expect_silent(
    formatted_CS_sim <- DAISIE:::DAISIE_format_CS(
    island_replicates = island_replicates,
    time = time,
    M = mainland_n,
    sample_freq = sample_freq,
    verbose = verbose
    )
  )
  expected_CS_format <- list()
  expected_CS_format[[1]] <- list()
  stt_all <- matrix(ncol = 5, nrow = 2)
  colnames(stt_all) <- c("Time", "nI", "nA", "nC", "present")
  stt_all[1, ] <- c(1, 0, 0, 0, 0)
  stt_all[2, ] <- c(0, 0, 0, 0, 0)
  expected_CS_format[[1]][[1]] <- list(island_age = 1,
                                  not_present = 1,
                                  stt_all = stt_all)
  expect_identical(formatted_CS_sim, expected_CS_format)
})

test_that("silent with non-empty island with correct output", {
  pars <- c(0.5, 0.1, 10, 1, 0.5)
  time <- 1
  mainland_n <- 1
  verbose <- FALSE
  sample_freq <- 1
  set.seed(1)
  island_replicates <- list()
  out <- list()
  out[[1]] <- DAISIE:::DAISIE_sim_core_cr(
    time = time,
    pars = pars,
    mainland_n = mainland_n,
    nonoceanic_pars = c(0, 0),
    hyper_pars = create_hyper_pars(d = 0, x = 0),
    area_pars = DAISIE::create_area_pars(
      max_area = 1,
      current_area = 1,
      proportional_peak_t = 0,
      total_island_age = 0,
      sea_level_amplitude = 0,
      sea_level_frequency = 0,
      island_gradient_angle = 0)
  )
  island_replicates[[1]] <- out
  expect_silent(
    formatted_CS_sim <- DAISIE:::DAISIE_format_CS( # nolint
      island_replicates = island_replicates,
      time = time,
      M = mainland_n,
      sample_freq = sample_freq,
      verbose = verbose
    )
  )
  expected_CS_format <- list()
  expected_CS_format[[1]] <- list()
  stt_all <- matrix(ncol = 5, nrow = 2)
  colnames(stt_all) <- c("Time", "nI", "nA", "nC", "present")
  stt_all[1, ] <- c(1, 0, 0, 0, 0)
  stt_all[2, ] <- c(0, 0, 0, 3, 1)
  expected_CS_format[[1]][[1]] <- list(island_age = 1,
                                       not_present = 0,
                                       stt_all = stt_all)
  expected_CS_format[[1]][[2]] <- list(branching_times = c(1.000000000000000,
                                                           0.244818166871655,
                                                           0.173128288990374,
                                                           0.029668240213840),
                                       stac = 2,
                                       missing_species = 0)
  expect_equal(formatted_CS_sim, expected_CS_format)
})

test_that("output with empty island and verbose = TRUE", {
  pars <- c(0, 1, 1, 0.0001, 0)
  time <- 1
  mainland_n <- 1
  verbose <- TRUE
  sample_freq <- 1
  set.seed(1)
  island_replicates <- list()
  out <- list()
  out[[1]] <- DAISIE:::DAISIE_sim_core_cr(
    time = time,
    pars = pars,
    mainland_n = mainland_n,
    nonoceanic_pars = c(0, 0),
    hyper_pars = create_hyper_pars(d = 0, x = 0),
    area_pars = DAISIE::create_area_pars(
      max_area = 1,
      current_area = 1,
      proportional_peak_t = 0,
      total_island_age = 0,
      sea_level_amplitude = 0,
      sea_level_frequency = 0,
      island_gradient_angle = 0)
  )
  island_replicates[[1]] <- out
  expect_message(
    formatted_CS_sim <- DAISIE:::DAISIE_format_CS(
      island_replicates = island_replicates,
      time = time,
      M = mainland_n,
      sample_freq = sample_freq,
      verbose = verbose
    ),
    regexp = "Island being formatted: 1/1"
  )
})

test_that("silent with non-empty 2 type island", {
  pars <- c(0.4, 0.1, 10, 1, 0.5, 0.4, 0.1, 10, 1, 0.5)
  total_time <- 1
  M <- 10
  mainland_n <- M
  verbose <- FALSE
  replicates <- 1
  sample_freq <- 10
  set.seed(1)
  island_replicates <- list()
  prop_type2_pool <- 0.4
  island_replicates <- DAISIE:::DAISIE_sim_min_type2(
    time = total_time,
    M = M,
    pars = pars,
    replicates = replicates,
    prop_type2_pool = prop_type2_pool,
    hyper_pars = create_hyper_pars(d = 0, x = 0),
    area_pars = DAISIE::create_area_pars(
      max_area = 1,
      current_area = 1,
      proportional_peak_t = 0,
      total_island_age = 0,
      sea_level_amplitude = 0,
      sea_level_frequency = 0,
      island_gradient_angle = 0),
    verbose = FALSE)
  expect_silent(
    formatted_CS_sim <- DAISIE:::DAISIE_format_CS(
      island_replicates = island_replicates,
      time = total_time,
      M = mainland_n,
      sample_freq = sample_freq,
      verbose = verbose
    )
  )
})

test_that("silent with non-empty 2 type island full stt", {
  pars <- c(0.4, 0.1, 10, 1, 0.5, 0.4, 0.1, 10, 1, 0.5)
  total_time <- 1
  M <- 10
  mainland_n <- M
  verbose <- FALSE
  replicates <- 1
  sample_freq <- Inf
  set.seed(1)
  island_replicates <- list()
  prop_type2_pool <- 0.4
  island_replicates <- DAISIE:::DAISIE_sim_min_type2(
    time = total_time,
    M = M,
    pars = pars,
    replicates = replicates,
    prop_type2_pool = prop_type2_pool,
    verbose = FALSE,
    hyper_pars = create_hyper_pars(d = 0, x = 0),
    area_pars = DAISIE::create_area_pars(
      max_area = 1,
      current_area = 1,
      proportional_peak_t = 0,
      total_island_age = 0,
      sea_level_amplitude = 0,
      sea_level_frequency = 0,
      island_gradient_angle = 0)
  )
  expect_silent(
    formatted_CS_sim <- DAISIE:::DAISIE_format_CS(
      island_replicates = island_replicates,
      time = total_time,
      M = mainland_n,
      verbose = verbose
    )
  )
})

test_that("abuse", {
  expect_error(
    DAISIE:::DAISIE_format_CS(
      island_replicates = "nonsense",
      time = "nonsense",
      M = "nonsense",
      sample_freq = "nonsense",
      verbose = "nonsense"
    )
  )
})


test_that("use full stt", {
  pars <- c(0.4, 0.2, 10, 2, 0.5)
  time <- 5
  mainland_n <- 1
  verbose <- FALSE
  sample_freq <- Inf
  set.seed(1)
  island_replicates <- list()
  out <- list()
  out[[1]] <- DAISIE:::DAISIE_sim_core_cr(
    time = time,
    pars = pars,
    mainland_n = mainland_n,
    nonoceanic_pars = c(0, 0),
    hyper_pars = create_hyper_pars(d = 0, x = 0),
    area_pars = DAISIE::create_area_pars(
      max_area = 1,
      current_area = 1,
      proportional_peak_t = 0,
      total_island_age = 0,
      sea_level_amplitude = 0,
      sea_level_frequency = 0,
      island_gradient_angle = 0)
  )
  island_replicates[[1]] <- out
  expect_silent(
    formatted_CS_sim <- DAISIE:::DAISIE_format_CS(
      island_replicates = island_replicates,
      time = time,
      M = mainland_n,
      sample_freq = sample_freq,
      verbose = verbose
    )
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$island_age,
    5
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$not_present,
    0
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$stt_all[2, ],
    c(Time = 4.62240908343582735, nI = 1.0, nA = 0.0, nC = 0.0, present = 1.0)
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$stt_all[5, ],
    c(Time = 3.81548257687248871, nI = 0.0, nA = 1.0, nC = 0.0, present = 1.0)
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$stt_all[11, ],
    c(Time = 2.22760715636035123, nI = 1.0, nA = 0.0, nC = 0.0, present = 1.0)
  )

  expect_equal(
    formatted_CS_sim[[1]][[2]]$branching_times,
    c(5.0000000000000000, 2.22760715636035, 1.3487418169725700, 0.0921013811906803)
  )

  expect_equal(
    formatted_CS_sim[[1]][[2]]$stac,
    3
  )

  expect_equal(
    formatted_CS_sim[[1]][[2]]$all_colonisations[[1]]$event_times,
    c(5.0, 2.22760716)
  )

  expect_equal(
    formatted_CS_sim[[1]][[2]]$all_colonisations[[1]]$species_type,
    "A"
  )
  expect_equal(
    formatted_CS_sim[[1]][[2]]$all_colonisations[[2]]$event_times,
    c(5.0, 1.34874181697257, 0.0921013811906803)
  )

  expect_equal(
    formatted_CS_sim[[1]][[2]]$all_colonisations[[2]]$species_type,
    "C"
  )

  expect_equal(
    formatted_CS_sim[[1]][[2]]$missing_species,
    0
  )
})

test_that("use complete stt with ontogeny", {
  total_time <- 10
  mainland_n <- 1
  verbose <- FALSE
  sample_freq <- Inf
  set.seed(2)
  island_replicates <- list()
  out <- list()
  pars = c(0.1, 0.5, 0.005, 1, 1)
  area_pars <- create_area_pars(
    max_area = 5000,
    current_area = 3000,
    proportional_peak_t = 0.5,
    total_island_age = 15,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0
  )
  hyper_pars <- create_hyper_pars(d = 0.2, x = 0.1)
  nonoceanic_pars <- c(0, 0)
  island_ontogeny = DAISIE:::translate_island_ontogeny("beta")
  sea_level = DAISIE:::translate_sea_level("const")
  peak <- DAISIE:::calc_peak(total_time = total_time,
                             area_pars = area_pars)
  Amax <- DAISIE:::get_global_max_area(total_time = total_time,
                                       area_pars = area_pars,
                                       peak = peak,
                                       island_ontogeny = island_ontogeny,
                                       sea_level = sea_level)
  Amin <- DAISIE:::get_global_min_area(total_time = total_time,
                                       area_pars = area_pars,
                                       peak = peak,
                                       island_ontogeny = island_ontogeny,
                                       sea_level = sea_level)
  out[[1]] <- DAISIE:::DAISIE_sim_core_time_dep(
    time = total_time,
    pars = pars,
    mainland_n = mainland_n,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    hyper_pars = hyper_pars,
    area_pars = area_pars,
    nonoceanic_pars = nonoceanic_pars,
    peak = peak,
    Amax = Amax,
    Amin = Amin,
    extcutoff = 1000
  )
  island_replicates[[1]] <- out
  expect_silent(
    formatted_CS_sim <- DAISIE:::DAISIE_format_CS(
      island_replicates = island_replicates,
      time = total_time,
      M = mainland_n,
      sample_freq = sample_freq,
      verbose = verbose
    )
  )

  expect_equal(
    formatted_CS_sim[[1]][[1]]$island_age,
    10
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$not_present,
    0
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$stt_all[5, ],
    c(Time = 6.2738223342025243, nI = 1.0, nA = 0.0, nC = 0.0, present = 1.0)
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$stt_all[12, ],
    c(Time = 3.2073998327641329, nI = 1.0, nA = 2.0, nC = 0.0, present = 1.0)
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$stt_all[16, ],
    c(Time = 2.1527893059827958, nI = 1.0, nA = 0.0, nC = 4.0, present = 1.0)
  )
  expect_equal(
    formatted_CS_sim[[1]][[2]]$branching_times,
    c(10.00000000000000, 6.27382233420251989, 4.99629097366935,
      2.14112644569665012, 1.13725415915635009, 0.54321135189076997,
      0.15055819094851100)
  )

  expect_equal(
    formatted_CS_sim[[1]][[2]]$stac,
    3
  )

  expect_equal(
    formatted_CS_sim[[1]][[2]]$missing_species,
    0
  )
})

test_that("full stt works with multiple replicates", {
  pars <- c(0.4, 0.2, 10, 2, 0.5)
  time <- 5
  mainland_n <- 1
  verbose <- FALSE
  sample_freq <- Inf
  set.seed(1)
  area_pars <- DAISIE::create_area_pars(
    max_area = 1,
    current_area = 1,
    proportional_peak_t = 0,
    total_island_age = 0,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0)
  hyper_pars <- create_hyper_pars(d = 0, x = 0)
  nonoceanic_pars <- c(0, 0)
  island_replicates <- list()
  island_replicates[[1]] <- island_replicates
  out <- list()
  out[[1]] <- DAISIE:::DAISIE_sim_core_cr(
    time = time,
    pars = pars,
    mainland_n = mainland_n,
    area_pars = area_pars,
    hyper_pars = hyper_pars,
    nonoceanic_pars = nonoceanic_pars
  )
  out[[2]] <- DAISIE:::DAISIE_sim_core_cr(
    time = time,
    pars = pars,
    mainland_n = mainland_n,
    area_pars = area_pars,
    hyper_pars = hyper_pars,
    nonoceanic_pars = nonoceanic_pars
  )
  island_replicates[[1]] <- out
  expect_silent(
    formatted_CS_sim <- DAISIE:::DAISIE_format_CS(
      island_replicates = island_replicates,
      time = time,
      M = mainland_n,
      sample_freq = sample_freq,
      verbose = verbose
    )
  )
})

test_that("full stt works with empty island", {
  pars <- c(0.4, 0.2, 10, 0.0000001, 0.5)
  total_time <- 1
  mainland_n <- 1
  verbose <- FALSE
  sample_freq <- Inf
  area_pars <- DAISIE::create_area_pars(
    max_area = 1,
    current_area = 1,
    proportional_peak_t = 0,
    total_island_age = 0,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0)
  hyper_pars <- create_hyper_pars(d = 0, x = 0)
  nonoceanic_pars <- c(0, 0)
  set.seed(1)
  replicates <- 2
  island_replicates <- list()
  for (rep in 1:replicates) {
    island_replicates[[rep]] <- list()
    full_list <- list()
    out <- list()
    for (m_spec in 1:mainland_n) {
      out$branching_times <- c(10)
        out <- DAISIE:::DAISIE_sim_core_cr(
          time = total_time,
          mainland_n = 1,
          pars = pars,
          area_pars = area_pars,
          hyper_pars = hyper_pars,
          nonoceanic_pars = nonoceanic_pars
        )
      full_list[[m_spec]] <- out
    }
    island_replicates[[rep]] <- full_list
  }
  expect_silent(
    formatted_CS_sim <- DAISIE:::DAISIE_format_CS(
      island_replicates = island_replicates,
      time = total_time,
      M = mainland_n,
      sample_freq = sample_freq,
      verbose = verbose
    )
  )
})

test_that("full stt with two trait states", {
  pars <- c(0.4, 0.2, 20, 2, 0.5)
  time <- 5
  mainland_n <- 0
  verbose <- FALSE
  sample_freq <- Inf
  trait_pars = create_trait_pars(
    trans_rate = 0,
    immig_rate2 = 1,
    ext_rate2 = 0.4,
    ana_rate2 = 0.8,
    clado_rate2 = 0.4,
    trans_rate2 = 0,
    M2 = 1)
  set.seed(1)
  island_replicates <- list()
  out <- list()
  out[[1]] <- DAISIE:::DAISIE_sim_core_trait_dependent(
    time = time,
    pars = pars,
    mainland_n = mainland_n,
    hyper_pars = create_hyper_pars(d = 0, x = 0),
    area_pars = DAISIE::create_area_pars(
      max_area = 1,
      current_area = 1,
      proportional_peak_t = 0,
      total_island_age = 0,
      sea_level_amplitude = 0,
      sea_level_frequency = 0,
      island_gradient_angle = 0),
    trait_pars = trait_pars
  )
  island_replicates[[1]] <- out
  expect_silent(
    formatted_CS_sim <- DAISIE:::DAISIE_format_CS(
      island_replicates = island_replicates,
      time = time,
      M = mainland_n,
      sample_freq = sample_freq,
      verbose = verbose,
      trait_pars = trait_pars
    )
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$island_age,
    5
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$not_present,
    0
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$stt_two_states[2, ],
    c(Time = 4.24481817, nI = 0.0, nA = 0.0, nC = 0.0, nI2 = 1.0, nA2 = 0.0, nC2 = 0.0, present = 1.0)
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$stt_two_states[5, ],
    c(Time = 3.61806444, nI = 0.0, nA = 0.0, nC = 0.0, nI2 = 1.0, nA2 = 0.0, nC2 = 2.0, present = 1.0)
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$stt_two_states[11, ],
    c(Time = 1.17170697, nI = 0.0, nA = 0.0, nC = 0.0, nI2 = 0.0, nA2 = 3.0, nC2 = 0.0, present = 1.0)
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$stt_all[2, ],
    c(Time = 4.24481817, nI = 1.0, nA = 0.0, nC = 0.0, present = 1.0)
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$stt_all[5, ],
    c(Time = 3.61806444, nI = 1.0, nA = 0.0, nC = 2.0, present = 1.0)
  )
  expect_equal(
    formatted_CS_sim[[1]][[1]]$stt_all[11, ],
    c(Time = 1.17170697, nI = 0.0, nA = 3.0, nC = 0.0, present = 1.0)
  )

  expect_equal(
    formatted_CS_sim[[1]][[2]]$branching_times,
    c(5, 4.24481817, 3.61806444, 1.25876988, 0.01277218)
  )

  expect_equal(
    formatted_CS_sim[[1]][[2]]$stac,
    3
  )

  expect_equal(
    formatted_CS_sim[[1]][[2]]$missing_species,
    0
  )
})

