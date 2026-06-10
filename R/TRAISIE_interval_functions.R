#' @keywords internal
dist_gamma_tma <- function(gamma,
                                   trait_mainland_ancestor,
                                   num_unique_states) {
  dist_gamma <- c()
  if (num_unique_states == 1) {

    dist_gamma <- gamma

  } else {
    if (sum(is.na(trait_mainland_ancestor)) ==
        length(trait_mainland_ancestor)) {
      dist_gamma <- gamma

    } else {
      s <- rep(0, num_unique_states)
      trait_index <- which(trait_mainland_ancestor == 1)
      s[trait_index] <- 1
      dist_gamma <- gamma * s
    }
  }
  return(dist_gamma)
}

###############
TRAISIE_use_stationary_weights <- function(Q) {

  diag(Q) <- 0
  diag(Q) <- -rowSums(Q)
  pi <- pracma::null(t(Q))
  diff <- 1
  for (i in 1:dim(pi)[[2]]) {
    if (pi[which.max(abs(pi[, i])), i] < 0) {
      pi[, i] <- -pi[, i]
    }
    if (any(pi[, i] < 0) && max(abs(pi[which(pi[, i] < 0), i])) > 1E-10) {
      warning("Substantial negative weights detected")
    }
    pi[which(pi[, i] < 0), i] <- 0
    pi[, i] <- pi[, i] / sum(pi[, i])
    diff_new <- abs(max(pi[, i])) - abs(min(pi[, i]))
    if (diff_new < diff) {
      diff <- diff_new
      i_choice <- i
    }
  }
  weight_states <- pi[, i_choice]
  return(weight_states)
}
################

#' @keywords internal
TRAISIE_compute_mainland_stationary_weights <- function(stat_weights,
                                                        Mp,
                                                        M,
                                                        num_hidden_states) {

  weights1 <- c()

  for (j in seq_along(Mp)) {

    idx <- ((j - 1) * num_hidden_states + 1):(j * num_hidden_states)

    weights_j <- stat_weights[idx]
    if (sum(weights_j) == 0) {
      weights_j <- weights_j
    } else {
      weights_j <- weights_j * (Mp[j] / M) / sum(weights_j)
    }
    weights1 <- c(weights1, weights_j)
  }

  weights1 <- weights1 / sum(weights1)

  weights2 <- stat_weights * (1 - (sum(Mp) / M)) / sum(stat_weights)

  weights <- weights1 + weights2
  weights <- weights / sum(weights)

  return(weights)
}

#' @keywords internal
TRAISIE_compute_likelihood_stationary_weights <- function(Lk_vec,
                                                          Mp,
                                                          M,
                                                          num_hidden_states) {

  weights1 <- c()

  for (j in seq_along(Mp)) {

    idx <- ((j - 1) * num_hidden_states + 1):(j * num_hidden_states)

    weights_j <- Lk_vec[idx]
    if (sum(weights_j) == 0) {

      weights_j <- weights_j
    } else {
      weights_j <- weights_j * (Mp[j] / M) / sum(weights_j)

    }
    weights1 <- c(weights1, weights_j)
  }

  weights1 <- weights1 / sum(weights1)

  weights2 <- Lk_vec * (1 - (sum(Mp) / M)) / sum(Lk_vec)

  weights <- weights1 + weights2
  weights <- weights / sum(weights)

  return(weights)
}

#' @keywords internal
TRAISIE_compute_mainland_weights <- function(Mp,
                                             M,
                                             num_hidden_states) {

  weights1 <- c()

  for (j in seq_along(Mp)) {

    weights_j <- rep((Mp[j] / M), num_hidden_states)
    weights1 <- c(weights1, weights_j)
  }

  weights2 <- (1 - (sum(Mp) / M))

  weights <- weights1 + weights2
  weights <- weights / sum(weights)

  return(weights)
}
#############



#' @keywords internal
TRAISIE_interval2 <- function(t, state, parameter) {
  with(as.list(c(state, parameter)), {
    lambdac <- parameter[[1]]
    mu      <- parameter[[2]]
    gamma   <- parameter[[3]]
    lambdaa <- parameter[[4]]
    q       <- parameter[[5]]
    p       <- parameter[[6]]
    trait_mainland_ancestor <- parameter[[7]]

    n <- (length(state) - 1) / 4

    dDE     <- numeric(n)
    dDM2    <- numeric(n)
    dDM3    <- numeric(n)
    dE      <- numeric(n)

    t_vec <- rowSums(q)

    DE  <- state[1:n]
    DM2 <- state[(n + 1):(n + n)]
    DM3 <- state[(n + n + 1):(n + n + n)]
    E   <- state[(n + n + n + 1):(n + n + n + n)]
    DA3 <- state[length(state)]

    q_mult_E   <- t(q %*% E)
    q_mult_DE  <- t(q %*% DE)
    q_mult_DM2 <- t(q %*% DM2)
    q_mult_DM3 <- t(q %*% DM3)


    # TODO: pass trait mainland_ancestor to interval functions
    # trait mainland ancestor is vector of probabilities

    dist_gamma <- dist_gamma_tma(gamma,
                                 trait_mainland_ancestor,
                                 n)

    dDE <- -(lambdac + mu + t_vec) * DE +
      2 * lambdac * DE * E +
      q_mult_DE

    dDM2 <- -(lambdac + mu + sum(dist_gamma) + lambdaa + t_vec) * DM2 +
      (lambdaa * DE + 2 * lambdac * DE * E + p * q_mult_DE) * DA3 +
      (1 - p) * q_mult_DM2

    dDM3 <- -(lambdac + mu + lambdaa + t_vec) * DM3 +
      (mu + lambdaa * E + lambdac * E * E + p * q_mult_E) * DA3 +
      (1 - p) * q_mult_DM3

    dE <- mu - (mu + lambdac + t_vec) * E +
      lambdac * E * E +
      q_mult_E

    dDA3 <- -sum(dist_gamma) * DA3 + sum(dist_gamma * DM3)
    return(list(c(dDE, dDM2, dDM3, dE, dDA3)))
  })
}

#' @keywords internal
TRAISIE_interval3 <- function(t, state, parameter) {
  with(as.list(c(state, parameter)), {
    lambdac <- parameter[[1]]
    mu      <- parameter[[2]]
    gamma   <- parameter[[3]]
    lambdaa <- parameter[[4]]
    q       <- parameter[[5]]
    p       <- parameter[[6]]
    trait_mainland_ancestor <- parameter[[7]]

    n <- (length(state) - 2) / 5

    dDE     <- numeric(n)
    dDM1    <- numeric(n)
    dDM2    <- numeric(n)
    dDM3    <- numeric(n)
    dE      <- numeric(n)

    t_vec <- rowSums(q)

    DE  <- state[1:n]
    DM1 <- state[(n + 1):(n + n)]
    DM2 <- state[(n + n + 1):(n + n + n)]
    DM3 <- state[(n + n + n + 1):(n + n + n + n)]
    E   <- state[(n + n + n + n + 1):(n + n + n + n + n)]
    DA2 <- state[length(state) - 1]
    DA3 <- state[length(state)]

    q_mult_E   <- t(q %*% E)
    q_mult_DE  <- t(q %*% DE)
    q_mult_DM1 <- t(q %*% DM1)
    q_mult_DM2 <- t(q %*% DM2)
    q_mult_DM3 <- t(q %*% DM3)

    dist_gamma <- dist_gamma_tma(gamma,
                                 trait_mainland_ancestor,
                                 n)


    dDE <- -(lambdac + mu + t_vec) * DE +
      2 * lambdac * DE * E +
      q_mult_DE

    dDM1 <- -(lambdac + mu + sum(dist_gamma) + lambdaa + t_vec) * DM1 +
      (mu + lambdaa * E + lambdac * E * E + p * q_mult_E) * DA2 +
      (1 - p) * q_mult_DM1 + sum(dist_gamma * DM2)

    dDM2 <- -(lambdac + mu + lambdaa + t_vec) * DM2 +
      (mu + lambdaa * E + lambdac * E * E + p * q_mult_E) * DA2 +
      (lambdaa * DE + 2 * lambdac * DE * E + p * q_mult_DE) * DA3 +
      (1 - p) * q_mult_DM2

    dDM3 <- -(lambdac + mu + lambdaa + t_vec) * DM3 +
      (mu + lambdaa * E + lambdac * E * E + p * q_mult_E) * DA3 +
      (1 - p) * q_mult_DM3

    dE <- mu - (mu + lambdac + t_vec) * E +
      lambdac * E * E +
      q_mult_E

    dDA2 <- -sum(dist_gamma) * DA2 + sum(dist_gamma * DM2)
    dDA3 <- -sum(dist_gamma) * DA3 + sum(dist_gamma * DM3)

    return(list(c(dDE, dDM1, dDM2, dDM3, dE, dDA2, dDA3)))
  })
}

#' @keywords internal
TRAISIE_interval4 <- function(t, state, parameter) {
  with(as.list(c(state, parameter)), {
    lambdac <- parameter[[1]]
    mu      <- parameter[[2]]
    gamma   <- parameter[[3]]
    lambdaa <- parameter[[4]]
    q       <- parameter[[5]]
    p       <- parameter[[6]]
    trait_mainland_ancestor <- parameter[[7]]

    n <- (length(state) - 1) / 2

    dDM1 <- numeric(n)
    dDE  <- numeric(n)

    t_vec <- rowSums(q)

    DM1 <- state[1:n]
    E   <- state[(n + 1):(n + n)]
    DA1 <- state[length(state)]

    q_mult_E   <- t(q %*% E)
    q_mult_DM1 <- t(q %*% DM1)


    dist_gamma <- dist_gamma_tma(gamma,
                                 trait_mainland_ancestor,
                                 n)

    dDM1 <- -(lambdac + mu + lambdaa + t_vec) * DM1 +
      (mu + lambdaa * E + lambdac * E * E + p * q_mult_E) * DA1 +
      (1 - p) * q_mult_DM1

    dE <- mu - (mu + lambdac + t_vec) * E +
      lambdac * E * E +
      q_mult_E

    dDA1 <- -sum(dist_gamma) * DA1 + sum(dist_gamma * DM1)

    return(list(c(dDM1, dE, dDA1)))
  })
}
