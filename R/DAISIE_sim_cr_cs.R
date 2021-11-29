#' Simulates island replicates with an clade-specific (CS) diversity-dependent
#' constant-rate process
#'
#' @inheritParams default_params_doc
#'
#' @return A list. The highest level of the least corresponds to each individual
#' replicate. See return for `DAISIE_sim_cr()` for details.
DAISIE_sim_cr_cs <- function(total_time,
                             M,
                             pars,
                             replicates,
                             nonoceanic_pars,
                             prop_type2_pool,
                             replicates_apply_type2,
                             sample_freq,
                             hyper_pars,
                             area_pars,
                             cond,
                             verbose) {
  island_replicates <- list()
  if (length(pars) == 5) {
    for (rep in 1:replicates) {
      island_replicates[[rep]] <- list()
      full_list <- list()
      if (cond == 0) {
        number_present <- -1
      } else {
        number_present <- 0
      }
      while (number_present < cond) {
        for (m_spec in 1:M) {
          full_list[[m_spec]] <- DAISIE_sim_core_cr(
            time = total_time,
            mainland_n = 1,
            pars = pars,
            nonoceanic_pars = nonoceanic_pars,
            hyper_pars = hyper_pars,
            area_pars = area_pars
          )
        }
        stac_vec <- unlist(full_list)[which(names(unlist(full_list)) == "stac")]
        present <- which(stac_vec != 0)
        number_present <- length(present)
      }
      island_replicates[[rep]] <- full_list
      if (verbose == TRUE) {
        print(paste("Island replicate ", rep, sep = ""))
      }
    }
  } else if (length(pars) == 10) {
    if (cond > 0) {
      warning(
        paste0(
          "Conditioning on number of colonisations is not implemented for 2
  type simulations. Returning result with no conditioning instead."
        )
      )
    }
    if (replicates_apply_type2 == TRUE) {
      island_replicates <- DAISIE_sim_min_type2(
        time = total_time,
        M = M,
        pars = pars,
        replicates = replicates,
        prop_type2_pool = prop_type2_pool,
        area_pars = area_pars,
        hyper_pars = hyper_pars,
        verbose = verbose)
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
          full_list[[m_spec]] <- DAISIE_sim_core_cr(
            time = total_time,
            mainland_n = 1,
            pars = c(lac_1,
                     mu_1,
                     K_1,
                     gam_1,
                     laa_1),
            area_pars = area_pars,
            hyper_pars = hyper_pars,
            nonoceanic_pars = c(0, 0))
          full_list[[m_spec]]$type1or2  <- 1
        }
        #### species of pool2
        for (m_spec in (pool1 + 1):(pool1 + pool2)) {
          full_list[[m_spec]] <- DAISIE_sim_core_cr(
            time = total_time,
            mainland_n = 1,
            pars = c(lac_2,
                     mu_2,
                     K_2,
                     gam_2,
                     laa_2),
            area_pars = area_pars,
            hyper_pars = hyper_pars,
            nonoceanic_pars = c(0, 0))
          full_list[[m_spec]]$type1or2 <- 2
        }
        island_replicates[[rep]] <- full_list
        if (verbose == TRUE) {
          print(paste("Island replicate ", rep, sep = ""))
        }
      }
    }
  }
  island_replicates <- DAISIE_format_CS(
    island_replicates = island_replicates,
    time = total_time,
    M = M,
    sample_freq = sample_freq,
    verbose = verbose
  )

  return(island_replicates)
}

