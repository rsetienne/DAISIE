DAISIE_sim_new <- function(
  time,
  M,
  pars,
  replicates,
  mainland_params = NULL,
  divdepmodel = "CS",
  ddmodel_sim = 11,
  island_type = "oceanic",
  nonoceanic_params = NULL,
  k_dist_params = NULL,
  num_guilds = NULL,
  prop_type2_pool = NA,
  replicates_apply_type2 = TRUE,
  sample_freq = 25,
  plot_sims = TRUE,
  island_ontogeny = "const",
  Apars = NULL,
  Epars = NULL,
  sea_level = "const",
  Spars = NULL,
  land_bridge = FALSE,
  land_bridge_threshold = NULL,
  keep_final_state = FALSE,
  stored_data = NULL,
  verbose = TRUE,
  ...
) {
  testit::assert(
    "island_ontogeny is not valid input. Specify 'const',\n
    'linear' or  ' beta'", is_island_ontogeny_input(island_ontogeny)
  )
  #TODO: TEST island_replicates INPUT! SANITIZE STORED_DATA INPUT! ASSERT + TEST
  if (!is.null(stored_data)) {
    start_midway <- TRUE
  } else {
    start_midway <- FALSE
  }
  # @richelbilderbeek
  if (!is.null(mainland_params)) {
    return(
      DAISIE_sim_with_mainland(
        time = time,
        M = M,
        pars = pars,
        replicates = replicates,
        mainland_params = mainland_params,
        divdepmodel = divdepmodel,
        prop_type2_pool = prop_type2_pool,
        replicates_apply_type2 = replicates_apply_type2,
        sample_freq = sample_freq
      )
    )
  }
  # Classic behavior
  totaltime <- time
  island_replicates <- list()
  if (divdepmodel == "IW") {
    if (length(pars) > 5) {
      stop("Island-wide carrying capacity model not yet implemented for
           two types of mainland species")
    }
    for (rep in 1:replicates) {
      island_replicates[[rep]] <- DAISIE_sim_core(
        time = totaltime,
        mainland_n = M,
        pars = pars,
        ddmodel_sim = ddmodel_sim,
        island_type = island_type,
        nonoceanic_params = nonoceanic_params,
        island_ontogeny = island_ontogeny,
        Apars = Apars,
        Epars = Epars,
        keep_final_state = keep_final_state,
        island_spec = NULL
      )
      if (verbose == TRUE) {
        print(paste("Island replicate ", rep, sep = ""))
      }
    }
    island_replicates <- DAISIE_format_IW(island_replicates = island_replicates,
                                          time = totaltime,
                                          M = M,
                                          sample_freq = sample_freq,
                                          island_type = island_type)
  }
  if (divdepmodel == "CS") {
    if (length(pars) == 5) {
      # Midway simulation
      if (!is.null(stored_data)) {
        for (rep in 1:replicates) {
          n_colonized_replicates <- length(stored_data[[rep]]) - 1
          colonized_island_spec <- list()
          for (k in 1:n_colonized_replicates) {
            colonized_island_spec[[k]] <- stored_data[[rep]][[k + 1]]$island_spec
          }
          island_replicates <- list()
          # Run each clade seperately
          full_list <- list()
          if (length(colonized_island_spec) > 0) {
            # Run midway clades
            for (m_spec in 1:n_colonized_replicates) {
              full_list[[m_spec]] <- DAISIE_sim_core(
                time = totaltime,
                mainland_n = 1,
                pars = pars,
                ddmodel_sim = ddmodel_sim,
                island_type = island_type,
                nonoceanic_params = nonoceanic_params,
                island_ontogeny = island_ontogeny,
                Apars = Apars,
                Epars = Epars,
                keep_final_state = keep_final_state,
                island_spec = colonized_island_spec[[m_spec]]
              )
            }
          } else {
            # Run empty clades that didn't get colonists
            for (m_spec in (n_colonized_replicates + 1):1000) {
              full_list[[m_spec]] <- DAISIE_sim_core(
                time = totaltime,
                mainland_n = 1,
                pars = pars,
                ddmodel_sim = ddmodel_sim,
                island_type = island_type,
                nonoceanic_params = nonoceanic_params,
                island_ontogeny = island_ontogeny,
                Apars = Apars,
                Epars = Epars,
                keep_final_state = keep_final_state,
                island_spec = NULL
              )
            }
          }
          island_replicates[[rep]] <- full_list
          if (verbose == TRUE) {
            print(paste("Island replicate ", rep, sep = ""))
          }
        }
      } else {
        # Simulation from empty island
        for (rep in 1:replicates) {
          island_replicates[[rep]] <- list()
          # Run each clade seperately
          full_list <- list()
          for (m_spec in 1:M) {
            full_list[[m_spec]] <- DAISIE_sim_core(
              time = totaltime,
              mainland_n = 1,
              pars = pars,
              ddmodel_sim = ddmodel_sim,
              island_type = island_type,
              nonoceanic_params = nonoceanic_params,
              k_dist_params = k_dist_params,
              island_ontogeny = island_ontogeny,
              Apars = Apars,
              Epars = Epars,
              keep_final_state = keep_final_state,
              island_spec = NULL
            )
          }
          island_replicates[[rep]] <- full_list
          if (verbose == TRUE) {
            print(paste("Island replicate ", rep, sep = ""))
          }
        }
      }
    }
    if (length(pars) == 10 & land_bridge == FALSE) {
      if (is.na(prop_type2_pool)) {
        stop("prop_type2_pool (fraction of mainland species that belongs to
             the second subset of species) must be specified when running
             model with two species types")
      }
      if (island_type == "nonoceanic") {
        stop("nonoceanic islands cannot have two type islands")
      }
      if (replicates_apply_type2 == TRUE) {
        island_replicates <- DAISIE_sim_min_type2(time = totaltime,
                                                  M = M,
                                                  pars = pars,
                                                  replicates = replicates,
                                                  prop_type2_pool = prop_type2_pool)
      } else {
        for (rep in 1:replicates) {
          pool2 <- DDD::roundn(M * prop_type2_pool)
          pool1 <- M - pool2
          lac_1 <- pars[1]
          mu_1 <- pars[2]
          K_1 <- pars[3]
          gam_1 <- pars[4]
          laa_1 <- pars[5]
          lac_2 <- pars[6]
          mu_2 <- pars[7]
          K_2 <- pars[8]
          gam_2 <- pars[9]
          laa_2 <- pars[10]
          full_list <- list()
          #### species of pool1
          for (m_spec in 1:pool1) {
            full_list[[m_spec]] <- DAISIE_sim_core(time = totaltime,
                                                   mainland_n = 1,
                                                   pars = c(lac_1,
                                                            mu_1,
                                                            K_1,
                                                            gam_1,
                                                            laa_1),
                                                   ddmodel_sim = ddmodel_sim)
            full_list[[m_spec]]$type1or2  <- 1
          }
          #### species of pool2
          for (m_spec in (pool1 + 1):(pool1 + pool2)) {
            full_list[[m_spec]] <- DAISIE_sim_core(time = totaltime,
                                                   mainland_n = 1,
                                                   pars = c(lac_2,
                                                            mu_2,
                                                            K_2,
                                                            gam_2,
                                                            laa_2),
                                                   ddmodel_sim = ddmodel_sim)
            full_list[[m_spec]]$type1or2 <- 2
          }
          island_replicates[[rep]] <- full_list
          if (verbose == TRUE) {
            print(paste("Island replicate ", rep, sep = ""))
          }
        }
      }
    }
    if (length(pars) == 10 & land_bridge == TRUE) {
      for (rep in 1:replicates) {
        island_replicates[[rep]] = list()
        full_list = list()
        parstmp <- c(pars[1:10], time - pars[11])
        for (m_spec in 1:M) {
          full_list[[m_spec]] = DAISIE_sim_core_lb(time = time, mainland_n = 1, parstmp)
        }
        island_replicates[[rep]] = full_list
        if (verbose == TRUE) {
          print(paste("Island replicate ", rep, sep = ""))
        }
      }
    }
    if (start_midway == TRUE) {
      island_replicates <- DAISIE_format_CS(
        island_replicates = island_replicates,
        time = totaltime,
        M = M,
        sample_freq = sample_freq,
        start_midway = start_midway,
        verbose = verbose)
    }
    island_replicates <- DAISIE_format_CS(island_replicates = island_replicates,
                                          time = totaltime,
                                          M = M,
                                          sample_freq = sample_freq,
                                          island_type = island_type,
                                          start_midway = start_midway,
                                          verbose = verbose)

  }

  if (divdepmodel == "GW") {
    if (!is.numeric(num_guilds)) {
      stop("num_guilds must be numeric")
    }
    if (length(pars) == 5) {
      guild_size <- M / num_guilds
      testit::assert(num_guilds < M)
      testit::assert(M %% num_guilds == 0)
      for (rep in 1:replicates) {
        island_replicates[[rep]] <- list()
        full_list <- list()
        for (m_spec in 1:num_guilds) {
          full_list[[m_spec]]  <- DAISIE_sim_core(time = totaltime,
                                                  mainland_n = guild_size,
                                                  pars = pars,
                                                  ddmodel_sim = ddmodel_sim,
                                                  island_type = island_type,
                                                  nonoceanic_params = nonoceanic_params,
                                                  k_dist_params = k_dist_params,
                                                  island_ontogeny = island_ontogeny,
                                                  Apars = Apars,
                                                  Epars = Epars,
                                                  keep_final_state = keep_final_state,
                                                  island_spec = NULL)
        }
        island_replicates[[rep]] <- full_list
        if (verbose == TRUE) {
          print(paste("Island replicate ", rep, sep = ""))
        }
      }
      island_replicates <- DAISIE_format_GW(island_replicates = island_replicates,
                                            time = totaltime,
                                            M = M,
                                            sample_freq = sample_freq,
                                            island_type = island_type,
                                            num_guilds = num_guilds,
                                            start_midway = start_midway,
                                            verbose = verbose)
    }
  }
  if (plot_sims == TRUE) {
    DAISIE_plot_sims(island_replicates)
  }
  return(island_replicates)
}
