get_ext_rate_dev <- function(mu,
                         hyper_pars,
                         extcutoff = 1000,
                         num_spec,
                         A,
                         trait_pars = NULL,
                         island_spec = NULL) {

  x <- hyper_pars$x
  if (is.null(trait_pars)) {
    ext_rate <- max(0, mu * (A ^ -x) * num_spec, na.rm = TRUE)
    ext_rate <- min(ext_rate, extcutoff, na.rm = TRUE)
    # testit::assert(ext_rate >= 0)
    return(ext_rate)
  } else {   ##species have two states
    if (is.matrix(island_spec) || is.null(island_spec)) {
      num_spec_trait1 <- length(which(island_spec[, 8] == "1"))
      num_spec_trait2 <- length(which(island_spec[, 8] == "2"))
    }
    ext_rate1 <- mu * num_spec_trait1
    ext_rate2 <- trait_pars$ext_rate2 * num_spec_trait2
    # testit::assert(is.numeric(ext_rate1))
    # testit::assert(is.numeric(ext_rate2))
    # testit::assert(ext_rate1 >= 0)
    # testit::assert(ext_rate2 >= 0)
    ext_list <- list(ext_rate1 = ext_rate1,
                     ext_rate2 = ext_rate2)
    return(ext_list)
  }
}

get_clado_rate_dev <- function(lac,
                           hyper_pars,
                           num_spec,
                           K,
                           A,
                           trait_pars = NULL,
                           island_spec = NULL) {
  # testit::assert(are_hyper_pars(hyper_pars))

  d <- hyper_pars$d
  if (is.null(trait_pars)) {
    clado_rate <- max(
      0, lac * num_spec * (A ^ d) * (1 - num_spec / (K * A)), na.rm = TRUE
    )
    # testit::assert(clado_rate >= 0)
    # testit::assert(is.numeric(clado_rate))
    return(clado_rate)
  }else{
    num_spec_trait1 <- length(which(island_spec[, 8] == "1"))
    num_spec_trait2 <- length(which(island_spec[, 8] == "2"))
    clado_rate1 <- max(
      0, lac * num_spec_trait1 * (1 - num_spec / K),
      na.rm = TRUE)
    clado_rate2 <- max(
      0, trait_pars$clado_rate2 * num_spec_trait2 * (1 - num_spec / K),
      na.rm = TRUE
    )
    # testit::assert(clado_rate1 >= 0)
    # testit::assert(clado_rate2 >= 0)
    # testit::assert(is.numeric(clado_rate1))
    # testit::assert(is.numeric(clado_rate2))
    clado_list <- list(clado_rate1 = clado_rate1,
                       clado_rate2 = clado_rate2)
    return(clado_list)
  }
}


get_immig_rate_dev <- function(gam,
                           A,
                           num_spec,
                           K,
                           mainland_n,
                           trait_pars = NULL,
                           island_spec = NULL) {

  if (is.null(trait_pars)) {
    immig_rate <- max(c(mainland_n * gam * (1 - (num_spec / (A * K))),
                        0), na.rm = TRUE)
    # testit::assert(is.numeric(immig_rate))
    # testit::assert(immig_rate >= 0)
    return(immig_rate)
  } else {
    mainland_n2 <- trait_pars$M2
    gam2 <- trait_pars$immig_rate2
    immig_rate1 <- max(c(mainland_n * gam * (1 - (num_spec / (A * K))),
                         0), na.rm = TRUE)
    immig_rate2 <- max(c(mainland_n2 * gam2 * (1 - (num_spec / (A * K))),
                         0), na.rm = TRUE)
    # testit::assert(is.numeric(immig_rate1))
    # testit::assert(immig_rate1 >= 0)
    # testit::assert(is.numeric(immig_rate2))
    # testit::assert(immig_rate2 >= 0)
    immig_list <- list(immig_rate1 = immig_rate1,
                       immig_rate2 = immig_rate2)
    return(immig_list)
  }
}
