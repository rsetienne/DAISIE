test_that("printing is done", {
  expect_message(
    print_parameters_and_loglik(
      pars = c(1, 2:6),
      loglik = -3,
      verbose = 3,
      type = 'clade_loglik'
    ),
    regexp = "Status of colonist: 1, Parameters: 2.000000, 3.000000, 4.000000, 5.000000, 6.000000, Loglikelihood: -3.000000")
  expect_message(print_parameters_and_loglik(pars = c(2:6), loglik = -3, verbose = 3, type = 'island_loglik'))

  expect_message(print_parameters_and_loglik(pars = c(2:6), loglik = -3, verbose = 3, type = 'island_ML'))
  expect_message(print_parameters_and_loglik(pars = c(1:11), loglik = D-3, verbose = 3, type = 'global_ML', distance_dep = 'power'))
  expect_message(print_parameters_and_loglik(pars = c(1:11), loglik = -3, verbose = 3, type = 'global_ML', distance_dep = 'sigmoidal_col'))
  expect_message(print_parameters_and_loglik(pars = data.frame(rbind(c(2:6),c(12:16))), loglik = -3, verbose = 3, type = 'multiple_island_ML'))

  expect_silent(print_parameters_and_loglik(pars = data.frame(rbind(c(2:6),c(12:16))), loglik = -3, verbose = 1, type = 'multiple_island_ML'))
  expect_silent(print_parameters_and_loglik(pars = data.frame(rbind(c(2:6),c(12:16))), loglik = -3, verbose = 0, type = 'multiple_island_ML'))
})
