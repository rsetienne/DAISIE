context("DAISIE_format_CS_trait")

test_that("sampled stt, 1 type, no geodynamics, oceanic, two trait states
          (same arguments as geodynamics, 5 pars)", {
            pars <- c(0.4, 0.2, 10, 2, 0.5)
            totaltime <- 1
            mainland_n <- 2
            verbose <- FALSE
            set.seed(1)
            replicates <- 3
            island_ontogeny = 0
            sea_level = 0
            extcutoff = 1000
            trait_pars <- create_trait_pars(
              trans_rate = 0,
              immig_rate2 = 2,
              ext_rate2 = 0.2,
              ana_rate2 = 0.5,
              clado_rate2 = 0.4,
              trans_rate2 = 0,
              M2 = 2)
            island_replicates <- list()
            verbose <- FALSE
            sample_freq <- 25

            for (rep in 1:replicates) {
              island_replicates[[rep]] <- list()
              full_list <- list()
              trait_pars_addcol <- create_trait_pars(trans_rate = 0,
                                                     immig_rate2 = 0,
                                                     ext_rate2 = 0,
                                                     ana_rate2 = 0,
                                                     clado_rate2 = 0,
                                                     trans_rate2 = 0,
                                                     M2 = 0)
              for (m_spec in 1:mainland_n) {
                full_list[[m_spec]] <- DAISIE_sim_core_trait_dependent(
                  time = totaltime,
                  mainland_n = 1,
                  pars = pars,
                  island_ontogeny = island_ontogeny,
                  sea_level = sea_level,
                  extcutoff = extcutoff,
                  hyper_pars = create_hyper_pars(d = 0, x = 0),
                  area_pars = DAISIE::create_area_pars(
                    max_area = 1,
                    current_area = 1,
                    proportional_peak_t = 0,
                    total_island_age = 0,
                    sea_level_amplitude = 0,
                    sea_level_frequency = 0,
                    island_gradient_angle = 0),
                  trait_pars = trait_pars_addcol
                )
              }
              for(m_spec in (mainland_n + 1):(mainland_n + trait_pars$M2))
              {
                trait_pars_onecolonize <- create_trait_pars(
                  trans_rate = trait_pars$trans_rate,
                  immig_rate2 = trait_pars$immig_rate2,
                  ext_rate2 = trait_pars$ext_rate2,
                  ana_rate2 = trait_pars$ana_rate2,
                  clado_rate2 = trait_pars$clado_rate2,
                  trans_rate2 = trait_pars$trans_rate2,
                  M2 = 1
                )
                full_list[[m_spec]] <- DAISIE_sim_core_trait_dependent(
                  time = totaltime,
                  mainland_n = 0,
                  pars = pars,
                  island_ontogeny = island_ontogeny,
                  sea_level = sea_level,
                  extcutoff = extcutoff,
                  hyper_pars = create_hyper_pars(d = 0, x = 0),
                  area_pars = DAISIE::create_area_pars(
                    max_area = 1,
                    current_area = 1,
                    proportional_peak_t = 0,
                    total_island_age = 0,
                    sea_level_amplitude = 0,
                    sea_level_frequency = 0,
                    island_gradient_angle = 0),
                  trait_pars = trait_pars_onecolonize
                )

              }
              island_replicates[[rep]] <- full_list
            }
            expect_silent(
              formatted_CS_sim <- DAISIE:::DAISIE_format_CS_sampled_stt(
                island_replicates = island_replicates,
                time = totaltime,
                M = mainland_n,
                verbose = verbose,
                sample_freq = sample_freq,
                trait_pars = trait_pars
              )
            )
          })
