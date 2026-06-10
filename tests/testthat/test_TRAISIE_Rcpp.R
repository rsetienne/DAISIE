test_that("R vs Rcpp", {
  phy <- ape::read.tree(text =
              "((A:2.1483291,B:2.1483291):0.209405,C:2.3577341):11.7431696;"
              )

  traits <- c(1, 1, 0)

  num_hidden_states <- 2
  num_observed_states <- 2
  num_total_states <- num_hidden_states * num_observed_states

  lambda_c <- rep(1, num_total_states)
  lambda_a <- rep(0, num_total_states)
  mus      <- rep(0, num_total_states)
  gammas   <- rep(0.01, num_total_states)
  p <- 0
  tma <- c(1, 0)

  parameters <- list()
  parameters[[1]] <- lambda_c
  parameters[[2]] <- mus
  parameters[[3]] <- gammas
  parameters[[4]] <- lambda_a
  parameters[[5]] <- NA    # placeholder
  parameters[[6]] <- p

  # we have to re-write the q rates as a matrix:
  base_matrix <- matrix(0, nrow = 2, ncol = 2)
  base_matrix[1, 2] <- 1
  base_matrix[2, 1] <- 1
  expanded_matrix <- secsse::q_doubletrans(traits = c(0, 1),
                                           masterBlock = base_matrix,
                                           diff.conceal = FALSE)

  expanded_matrix[expanded_matrix == 1] <- 0.1


  parameters[[5]] <- expanded_matrix

  all_tma <- list(c(1, 1, 0, 0),
                  c(0, 0, 1, 1),
                  c(NA, NA, NA, NA))

  testthat::expect_warning(
  for (i in seq_along(all_tma)) {
    # this call generates warnings, but they are harmless
    res_hidden <- TRAISIE::loglik_R_tree(parameter = parameters,
                                        phy = phy,
                                        traits = traits,
                                        trait_mainland_ancestor = all_tma[[i]],
                                        num_hidden_states = 2,
                                        sampling_fraction = c(1, 1))

    res_cpp <- TRAISIE::loglik_cpp_tree(parameter = parameters,
                                       phy = phy,
                                       traits = traits,
                                       num_hidden_states = 2,
                                       sampling_fraction = c(1, 1),
                                       trait_mainland_ancestor = all_tma[[i]])
    testthat::expect_equal(res_hidden, res_cpp)
  }
  )
})
