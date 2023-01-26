test_that("printing is done", {
  print_parameters_and_loglik(pars = c(1,2:6), loglik = -3, verbose = T, type = 'clade_loglik')
  print_parameters_and_loglik(pars = c(2:6), loglik = -3, verbose = T, type = 'island_loglik')
  print_parameters_and_loglik(pars = c(2:6), loglik = -3, verbose = T, type = 'island_ML')
  print_parameters_and_loglik(pars = c(1:11), loglik = -3, verbose = T, type = 'global_ML', distance_dep = 'power')
  print_parameters_and_loglik(pars = c(1:11), loglik = -3, verbose = T, type = 'global_ML', distance_dep = 'sigmoidal_col')
  print_parameters_and_loglik(pars = data.frame(rbind(c(2:6),c(12:16))), loglik = -3, verbose = T, type = 'multiple_island_ML')
})
