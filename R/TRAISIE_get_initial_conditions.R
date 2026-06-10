#' @keywords internal
TRAISIE_initial_cond2_one_state <- function(status,
                                    res,
                                    trait,
                                    num_observed_states,
                                    num_hidden_states,
                                    brts,
                                    sampling_fraction,
                                    trait_mainland_ancestor) {

  initial_conditions2 <- c()

  DE  <- rep(0, 1)
  DM2 <- rep(0, 1)
  DM3 <- rep(0, 1)
  E   <- rep(0, 1)
  DA3 <- 1
  if (status == 2 && length(brts) > 2 || status == 3 && length(brts) > 2) {
    initial_conditions2 <- c(res[1],                      ## DE
                            (res[1]) * res[length(res)], ## DM2
                             res[2],          ## DM3
                             res[3],  ## E
                             res[length(res)])              ## DA3
    # pre-emptive return because this one is constructed differently.
    return(matrix(initial_conditions2, nrow = 1))
  } else if (status == 2 && length(brts) == 2) {

    DE[1] <- sampling_fraction[1]
    E[1]  <- 1 - sampling_fraction[1]

  } else if (status == 3 && length(brts) == 2) {
    DE[1] <- sampling_fraction[1]
    E[1] <- 1 - sampling_fraction[1]
    DM3[1] <- 1
    DA3 <- 0
  } else if (status == 4) {

    DM2[1] <- 1

  } else if (status == 8) {

    DM2[1] <- 1

  } else if (status == 9)  {


    DE[1] <- sampling_fraction[1]
    E[1]  <- 1 - sampling_fraction[1]

  }
  initial_conditions2 <- c(DE, DM2, DM3, E, DA3)
  return(matrix(initial_conditions2, nrow = 1))
}


# initial conditions system of equation interval 2
#' @keywords internal
TRAISIE_get_initial_conditions2 <- function(status,
                                    res,
                                    trait,
                                    num_observed_states,
                                    num_hidden_states,
                                    brts,
                                    sampling_fraction,
                                    trait_mainland_ancestor) {
  n <- num_observed_states * num_hidden_states
  num_unique_states <- n

  initial_conditions2 <- c()

  DE  <- rep(0, num_unique_states)
  DM2 <- rep(0, num_unique_states)
  DM3 <- rep(0, num_unique_states)
  E   <- rep(0, num_unique_states)
  DA3 <- 1

  # only use the sampling fraction of the focal trait, assuming traits start
  # counting at 0.


  if (num_unique_states == 1) {
    return(initial_cond2_one_state(status,
                                   res,
                                   trait,
                                   num_observed_states,
                                   num_hidden_states,
                                   brts,
                                   sampling_fraction,
                                   trait_mainland_ancestor))
  }  else {
    if (status == 2 && length(brts) > 2 || status == 3 && length(brts) > 2) {
      initial_conditions2 <- c(res[1:n],                      ## DE
                              (res[1:n]) * res[length(res)], ## DM2
                               res[(n + 1):(n + n)],          ## DM3
                               res[(n + n + 1):(n + n + n)],  ## E
                               res[length(res)])              ## DA3
      # pre-emptive return because this one is constructed differently.
      return(matrix(initial_conditions2, nrow = 1))
    } else if (status == 2 && length(brts) == 2) {

      # secret assumption of trait being a single value
      if (length(trait) > 1) {
        stop("status == 2 assumes trait to be single value, found vector")
      }
      if (is.na(trait)) {

        s <- c()
        for (i in seq_along(sampling_fraction)) {
          s <- c(s, rep(sampling_fraction[i], num_hidden_states))
        }

        DE[1:n] <- s
        E[1:n]  <- 1 - s

        # Only apply the change if NOT all s are 1
        if (!all(s == 1)) {
          DE[1:n][s == 1] <- 0
        }
      } else {


        DE[(num_hidden_states * trait + 1):
             (num_hidden_states + trait * num_hidden_states)] <- sampling_fraction[1 + trait]
        E[(num_hidden_states * trait + 1):
             (num_hidden_states + trait * num_hidden_states)] <- 1 - sampling_fraction[1 + trait]

        rest_idx <- setdiff(seq_along(E), (num_hidden_states * trait + 1):(num_hidden_states + num_hidden_states * trait))
        for (i in rest_idx) {
          trait_i <- (i - 1) %/% num_hidden_states
          sf_i <- sampling_fraction[1 + trait_i]
          E[i] <- if (sf_i == 1) 0 else 1 - sf_i
        }
      }

    } else if (status == 3 && length(brts) == 2) {

     DE[(num_hidden_states * trait + 1):
           (num_hidden_states + trait * num_hidden_states)] <- sampling_fraction[1 + trait]
      E[(num_hidden_states * trait + 1):
          (num_hidden_states + trait * num_hidden_states)] <- 1 - sampling_fraction[1 + trait]


      DM3[(num_hidden_states * trait_mainland_ancestor + 1):
            (num_hidden_states + trait_mainland_ancestor * num_hidden_states)] <- 1
        rest_idx <- setdiff(seq_along(E), (num_hidden_states * trait + 1):(num_hidden_states + num_hidden_states * trait))
        for (i in rest_idx) {
          trait_i <- (i - 1) %/% num_hidden_states
          sf_i <- sampling_fraction[1 + trait_i]
          E[i] <- ifelse(sf_i == 1, 0, 1 - sf_i)
        }
    } else if (status == 4) {

      if (length(trait) > 1) {
        stop("status == 4 assumes trait to be single value, found vector")
      }

      if (is.na(trait)) {

        DM2[1:n] <- 1

      } else {
        DM2[(num_hidden_states * trait + 1):
              (num_hidden_states + trait * num_hidden_states)] <- 1
      }

    } else if (status == 8) {
      if (length(trait) > 1) {
        stop("status == 8 assumes trait to be single value, found vector")
      }

      if (is.na(trait)) {

        DM2[1:n] <- 1

      } else {
        DM2[(num_hidden_states * trait + 1):
              (num_hidden_states + trait * num_hidden_states)] <- 1
      }
    } else if (status == 9)  {
      if (length(trait) > 1) {
        stop("status == 9 assumes trait to be single value, found vector")
      }

      if (is.na(trait)) {

        s <- c()
        for (i in seq_along(sampling_fraction)) {
          s <- c(s, rep(sampling_fraction[i], num_hidden_states))
        }

        DE[1:n] <- s
        E[1:n]  <- 1 - s

        # Only apply the change if NOT all s are 1
        if (!all(s == 1)) {
          DE[1:n][s == 1] <- 0
        }
      }else  {

        DE[(num_hidden_states * trait + 1):
             (num_hidden_states + trait * num_hidden_states)] <- sampling_fraction[1 + trait]
        E[(num_hidden_states * trait + 1):
            (num_hidden_states + trait * num_hidden_states)] <- 1 - sampling_fraction[1 + trait]

        rest_idx <- setdiff(seq_along(E), (num_hidden_states * trait + 1):(num_hidden_states + num_hidden_states * trait))
        for (i in rest_idx) {
          trait_i <- (i - 1) %/% num_hidden_states
          sf_i <- sampling_fraction[1 + trait_i]
          E[i] <- if (sf_i == 1) 0 else 1 - sf_i
        }
      }
    }
  }
  initial_conditions2 <- c(DE, DM2, DM3, E, DA3)
  return(matrix(initial_conditions2, nrow = 1))
}


# initial conditions system of equation interval 3
#' @keywords internal
TRAISIE_get_initial_conditions3 <- function(status,
                                    res,
                                    num_observed_states,
                                    num_hidden_states,
                                    trait,
                                    sampling_fraction,
                                    solution) {

  n <- num_observed_states * num_hidden_states
  num_unique_states <- n

  initial_conditions3 <- c()
  DE  <- rep(0, num_unique_states)
  DM2 <- rep(0, num_unique_states)
  DM1 <- rep(0, num_unique_states)
  DM3 <- rep(0, num_unique_states)
  E   <- rep(0, num_unique_states)
  DA2 <- 0
  DA3 <- 1


  # only use the sampling fraction of the focal trait, assuming traits start
  # counting at 0.

  if (num_unique_states == 1) {


    if (status == 1) {


      DM2[1:n] <- sampling_fraction


      initial_conditions3 <- c(DE, DM1, DM2, DM3, E, DA2, DA3)
    } else if (status == 6) {
      initial_conditions3 <- c(res[1:n],                                              ## DE
                               rep(0, n),                                             ## DM1
                               (res[1:n]) * res[length(res)],                         ## DM2
                               res[(n + 1):(n + n)],                                  ## DM3
                               res[(n + n + 1):(n + n + n)],                          ## E
                               0,                                                     ## DA2
                               res[length(res)])                                      ## DA3
    } else if (status == 5) {

      DE[1:n] <- sampling_fraction[1:n]
      E[1:n]  <- 1 - sampling_fraction[1:n]


      initial_conditions3 <- c(DE, DM1, DM2, DM3, E, DA2, DA3)
    } else if (status == 8 || status == 9) {
      initial_conditions3 <- c(solution[2, ][1:n],
                               rep(0, n),                                             ### DE: select DE in solution2
                               solution[2, ][(n + 1):(n + n)],                        ### DM2: select DM2 in solution2
                               solution[2, ][(n + n + 1):(n + n + n)],                ### DM3: select DM3 in solution2
                               solution[2, ][(n + n + n + 1):(n + n + n + n)],        ### E: select E in solution2
                               0,
                               solution[2, ][length(solution[2, ])])                 ### DA3: select DA3 in solution2
    }

  } else {
    if (status == 1) {
      if (length(trait) > 1) {
        stop("status == 1 assumes trait to be single value, found vector")
      }


      if (is.na(trait)) {

        DM2[1:n] <- 1

      } else if (trait == trait) {
        DM2[(num_hidden_states * trait + 1):
              (num_hidden_states + trait * num_hidden_states)] <- sampling_fraction[1 + trait]
      }

      initial_conditions3 <- c(DE, DM1, DM2, DM3, E, DA2, DA3)
    } else if (status == 6) {
      initial_conditions3 <- c(res[1:n],                                              ## DE
                               rep(0, n),                                             ## DM1
                               (res[1:n]) * res[length(res)],                         ## DM2
                               res[(n + 1):(n + n)],                                  ## DM3
                               res[(n + n + 1):(n + n + n)],                          ## E
                               0,                                                     ## DA2
                               res[length(res)])                                      ## DA3
    } else if (status == 5) {
      if (length(trait) > 1) {
        stop("status == 5 assumes trait to be single value, found vector")
      }

      if (is.na(trait)) {

        s <- c()
        for (i in seq_along(sampling_fraction)) {
          s <- c(s, rep(sampling_fraction[i], num_hidden_states))
        }

        DE[1:n] <- s
        E[1:n]  <- 1 - s

        # Only apply the change if NOT all s are 1
        if (!all(s == 1)) {
          DE[1:n][s == 1] <- 0
        }
      }else  {
        DE[(num_hidden_states * trait + 1):
             (num_hidden_states + trait * num_hidden_states)] <- sampling_fraction[1 + trait]
        E[(num_hidden_states * trait + 1):
            (num_hidden_states + trait * num_hidden_states)] <- 1 - sampling_fraction[1 + trait]

        rest_idx <- setdiff(seq_along(E), (num_hidden_states * trait + 1):(num_hidden_states + num_hidden_states * trait))
        for (i in rest_idx) {
          trait_i <- (i - 1) %/% num_hidden_states
          sf_i <- sampling_fraction[1 + trait_i]
          E[i] <- if (sf_i == 1) 0 else 1 - sf_i
        }
      }

      initial_conditions3 <- c(DE, DM1, DM2, DM3, E, DA2, DA3)
    } else if (status == 8 || status == 9) {
      initial_conditions3 <- c(solution[2, ][1:n],
                               rep(0, n),                                             ### DE: select DE in solution2
                               solution[2, ][(n + 1):(n + n)],                        ### DM2: select DM2 in solution2
                               solution[2, ][(n + n + 1):(n + n + n)],                ### DM3: select DM3 in solution2
                               solution[2, ][(n + n + n + 1):(n + n + n + n)],        ### E: select E in solution2
                               0,
                               solution[2, ][length(solution[2, ])])                 ### DA3: select DA3 in solution2
    }
  }
  return(matrix(initial_conditions3, nrow = 1))
}


# initial conditions system of equation interval 4
#' @keywords internal
TRAISIE_get_initial_conditions4 <- function(status,
                                    solution,
                                    parameter,
                                    trait_mainland_ancestor,
                                    num_observed_states,
                                    num_hidden_states) {

  n <- num_observed_states * num_hidden_states
  num_unique_states <- n

  if (num_unique_states == 1) {

    if (status == 2 || status == 3 || status == 4) {

      dist_gamma <- dist_gamma_tma(parameter[[3]],
                                   trait_mainland_ancestor,
                                   n)

      initial_conditions4 <- c(rep(sum(dist_gamma * (solution[2, ][(n + 1):(n + n)])), n), ### DM1: select DM2 in solution2
                               solution[2, ][(n + n + n + 1):(n + n + n + n)],                                                                       ### E: select E in solution2
                               sum(dist_gamma * (solution[2, ][(n + 1):(n + n)])))          ### DA1: select DM2 in solution2

    } else if (status == 1 || status == 5 || status == 6) {
      initial_conditions4 <- c(solution[2, ][(n + 1):(n + n)],                                 ### DM1: select DM1 in solution1
                               solution[2, ][(n + n + n + n + 1):(n + n + n + n + n)],         ### E: select E in solution1
                               solution[2, ][length(solution[2, ]) - 1])                        ### DA1: select DA2 in solution1

    } else if (status == 8 || status == 9) {
      initial_conditions4 <- c(solution[2, ][(n + 1):(n + n)],                                 ### DM1: select DM2 in solution3
                               solution[2, ][(n + n + n + n + 1):(n + n + n + n + n)],         ### E: select E in solution3
                               solution[2, ][length(solution[2, ]) - 1])                       ### DA1: select DA2 in solution3

    }


  } else {
    if (status == 2 || status == 3 || status == 4) {
      if (all(is.na(trait_mainland_ancestor))) {

        dist_gamma <- dist_gamma_tma(parameter[[3]],
                                     trait_mainland_ancestor,
                                     n)
        initial_conditions4 <- c(rep(sum(dist_gamma * (solution[2, ][(n + 1):(n + n)])), n), ### DM1: select DM2 in solution2
                                 solution[2, ][(n + n + n + 1):(n + n + n + n)],         ### E: select E in solution2
                                 sum(dist_gamma * (solution[2, ][(n + 1):(n + n)])))          ### DA1: select DM2 in solution2
      } else {
        #if the trait state of the species at the stem is known

        dist_gamma <- dist_gamma_tma(parameter[[3]],
                                     trait_mainland_ancestor,
                                     n)

        initial_conditions4 <- c(rep(sum(dist_gamma * (solution[2, ][(n + 1):(n + n)])), n), ### DM1: select DM2 in solution2
                                 solution[2, ][(n + n + n + 1):(n + n + n + n)],                                                                       ### E: select E in solution2
                                 sum(dist_gamma * (solution[2, ][(n + 1):(n + n)])))          ### DA1: select DM2 in solution2
      }
    } else if (status == 1 || status == 5 || status == 6) {
      initial_conditions4 <- c(solution[2, ][(n + 1):(n + n)],                                 ### DM1: select DM1 in solution1
                               solution[2, ][(n + n + n + n + 1):(n + n + n + n + n)],         ### E: select E in solution1
                               solution[2, ][length(solution[2, ]) - 1])                        ### DA1: select DA2 in solution1

    } else if (status == 8 || status == 9) {
      initial_conditions4 <- c(solution[2, ][(n + 1):(n + n)],                                 ### DM1: select DM2 in solution3
                               solution[2, ][(n + n + n + n + 1):(n + n + n + n + n)],         ### E: select E in solution3
                               solution[2, ][length(solution[2, ]) - 1])                       ### DA1: select DA2 in solution3

    }

  }
  return(matrix(initial_conditions4, nrow = 1))
}
