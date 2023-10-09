test_that("print_parameters_and_loglik works", {
  testthat::expect_message(
    print_parameters_and_loglik(
      pars = c(1, 2:6),
      loglik = -3,
      verbose = 3,
      type = 'clade_loglik'
    ),
    regexp = "Status of colonist: 1, Parameters: 2.000000, 3.000000, 4.000000, 5.000000, 6.000000, Loglikelihood: -3.000000")
  testthat::expect_message(print_parameters_and_loglik(pars = c(2:6), loglik = -3, verbose = 3, type = 'island_loglik'))

  testthat::expect_message(print_parameters_and_loglik(pars = c(2:6), loglik = -3, verbose = 3, type = 'island_ML'))
  testthat::expect_message(print_parameters_and_loglik(pars = c(1:11), loglik = -3, verbose = 3, type = 'global_ML', distance_dep = 'power'))
  testthat::expect_message(print_parameters_and_loglik(pars = c(1:11), loglik = -3, verbose = 3, type = 'global_ML', distance_dep = 'sigmoidal_col'))
  testthat::expect_message(print_parameters_and_loglik(pars = data.frame(rbind(c(2:6),c(12:16))), loglik = -3, verbose = 3, type = 'multiple_island_ML'))

  testthat::expect_silent(print_parameters_and_loglik(pars = data.frame(rbind(c(2:6),c(12:16))), loglik = -3, verbose = 1, type = 'multiple_island_ML'))
  testthat::expect_silent(print_parameters_and_loglik(pars = data.frame(rbind(c(2:6),c(12:16))), loglik = -3, verbose = 0, type = 'multiple_island_ML'))
})

test_that("print_ml_par_settings works", {

  testthat::expect_message(
    print_ml_par_settings(
      namepars = letters[1:5],
      idparsopt = 1:5,
      idparsfix = NULL,
      idparsnoshift = 1:5,
      all_no_shift = 1:5,
      verbose = 1
    ),
    regexp = "You are optimizing: a b c d e
You are fixing: nothing")

  testthat::expect_message(
    print_ml_par_settings(
      namepars = letters[1:5],
      idparsopt = 1:4,
      idparsfix = 5,
      idparsnoshift = 1:5,
      all_no_shift = 1:5,
      verbose = 1
    ),
    regexp = "You are optimizing: a b c d
You are fixing: e")

  testthat::expect_message(
    print_ml_par_settings(
      namepars = letters[1:5],
      idparsopt = 1:5,
      idparsfix = NULL,
      idparsnoshift = 1:4,
      all_no_shift = 1:4,
      verbose = 1
    ),
    regexp = "You are optimizing: a b c d e
You are fixing: nothing
You are not shifting: a b c d")

  testthat::expect_silent(
    print_ml_par_settings(
      namepars = letters[1:5],
      idparsopt = 1:5,
      idparsfix = NULL,
      idparsnoshift = 1:4,
      all_no_shift = 1:4,
      verbose = 0
    )
  )

})

test_that("print_init_ll works", {
  testthat::expect_message(
    print_init_ll(
      initloglik = 10,
      verbose = 1
    ),
    regexp = "Calculating the likelihood for the initial parameters.
The loglikelihood for the initial parameter values is 10
Optimizing the likelihood - this may take a while.")

  testthat::expect_silent(
    print_init_ll(
      initloglik = 10,
      verbose = 0
    )
  )
})
