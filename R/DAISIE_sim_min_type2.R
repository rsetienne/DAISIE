#' Internal function of the DAISIE simulation
#'
#' @inheritParams default_params_doc
#' @keywords internal
DAISIE_sim_min_type2 <- function(time,
                                 M,
                                 pars,
                                 replicates,
                                 prop_type2_pool,
                                 area_pars,
                                 hyper_pars,
                                 verbose = TRUE) {
  island_replicates <- list()
  n_islands_with_type2 <- 0
  counter <- 0
  while (n_islands_with_type2 < replicates) {
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
      full_list[[m_spec]] <- DAISIE_sim_core_cr(
        time = time,
        mainland_n = 1,
        pars = c(lac_1,
                 mu_1,
                 K_1,
                 gam_1,
                 laa_1
        ),
        area_pars = area_pars,
        hyper_pars = hyper_pars,
        nonoceanic_pars = c(0, 0)
      )
      full_list[[m_spec]]$type1or2 <- 1
    }
    #### species of pool2
    for (m_spec in (pool1 + 1):(pool1 + pool2)) {
      full_list[[m_spec]] <- DAISIE_sim_core_cr(
        time = time,
        mainland_n = 1,
        pars = c(lac_2, mu_2, K_2, gam_2, laa_2),
        area_pars = area_pars,
        hyper_pars = hyper_pars,
        nonoceanic_pars = c(0, 0))
      full_list[[m_spec]]$type1or2 <- 2
    }
    type_2s <- which(unlist(full_list)[which(names(unlist(full_list)) == "type1or2")] == 2)
    number_type2_species_colonized <- 0
    for (i in 1:length(type_2s)) {
      the_row <- full_list[[type_2s[i]]]$stt_table[nrow(full_list[[type_2s[i]]]$stt_table), ]
      if (sum(the_row) > 0) {
        number_type2_species_colonized <- number_type2_species_colonized + 1
      }
    }
    if (number_type2_species_colonized > 0) {
      n_islands_with_type2 <- n_islands_with_type2 + 1
      island_replicates[[length(island_replicates) + 1]] <- list()
      island_replicates[[length(island_replicates)]] <- full_list
      if (verbose == TRUE) {
        message(
          "Number of island replicates with type 2 species: ",
          length(island_replicates)
        )
      }
    }
    counter <- counter + 1
    if (verbose == TRUE) {
      message("Island ", counter)
    }
  }
  return(island_replicates)
}
