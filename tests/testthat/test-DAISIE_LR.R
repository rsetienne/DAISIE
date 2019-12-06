test_that("bootstrapping is silent", {
  skip("WIP")
  utils::data(Galapagos_datalist)
  expect_silent(DAISIE_LR(datalist = Galapagos_datalist,
                          initparsopt_dd = c(2.5, 2.7, 20, 0.009, 1.01),
                          initparsopt_di = c(2.5, 2.7, 0.009, 1.01)))
})

test_that("calculates maximum likelihood for data for each model", {
  skip("WIP")
  expect_true()
})


# context("use step 1 of bootstrapping")
#
# test_that("run ML estimation under DI and DD and calculate likelihood ratio", {
#   utils::data(Galapagos_datalist)
#
#   #not sure these tests work
#     expect_silent(DAISIE_LR(
#       datalist = Galapagos_datalist,
#       initparsopt_di = c(2.5, 2.7, 0.009, 1.01),
#       initparsopt_dd = c(2.5,2.7,20,0.009,1.01),
#       idparsopt = 1:5,
#       parsfix = NULL,
#       idparsfix = NULL,
#       ddmodel = 11))
#     #not sure these tests work
#     expect_message(DAISIE_LR(
#       datalist = Galapagos_datalist,
#       initparsopt_di = c(2.5, 2.7, 0.009, 1.01),
#       initparsopt_dd = c(2.5,2.7,20,0.009,1.01),
#       idparsopt = 1:5,
#       parsfix = NULL,
#       idparsfix = NULL,
#       ddmodel = 11),
#       c("Estimating parameters under the diversity-independent model \n",
#       "Estimating parameters under the diversity-dependent model \n"),
#       all = TRUE
#       )
# })
#
#
# context("use step 2 of bootstrapping")
#
# test_that("run DI simulation from DI MLE", {
#
#   output <- DAISIE_LR(
#     datalist = Galapagos_datalist,
#     initparsopt_di = c(2.5, 2.7, 0.009, 1.01),
#     initparsopt_dd = c(2.5,2.7,20,0.009,1.01),
#     idparsopt = 1:5,
#     parsfix = NULL,
#     idparsfix = NULL,
#     ddmodel = 11,
#     endmc = 1)
#   expect_true(is.list(output)
# }
#
# context("use step 3 of bootstrapping")
#
# test_that("DI ML estimation for DI and DD simulations", {
#
#   output <- DAISIE_LR(
#     datalist = Galapagos_datalist,
#     initparsopt_di = c(2.5, 2.7, 0.009, 1.01),
#     initparsopt_dd = c(2.5,2.7,20,0.009,1.01),
#     idparsopt = 1:5,
#     parsfix = NULL,
#     idparsfix = NULL,
#     ddmodel = 11,
#     endmc = 5)
#   expect_true(is.list(output)
#
# }
#
# #Step 4
# likelihood_ratio_zero
# likelihood_ratio[mc]
#
# pvalue <- calc_p_value(out$LR[2:(endmc + 1)],out$LR[1])
#
# calc_p_value <- function(samples, x) {
#   samplessort = sort(samples)
#   pup = which(samplessort > x)
#   p_value = (length(pup) + 1)/ (length(samples) + 1)
#   return(p_value)
# }
#
# #Step 5
# likelihood_ratio_alpha <- as.numeric(stats::quantile(out$LR[2:(endmc + 1)],
#                                                      1 - alpha,
#                                                      type = 4))
#
# #Step 6
# dd_sims <- list()
# time <- datalist[[1]]$island_age
# M <- datalist[[1]]$not_present + (length(datalist) - 1)
# init_dd_pars <- as.numeric(c(init_dd_ml[1:5]))
#
#
# cat("\nSimulating under DD \n")
# for (me in 1:endmc) {
#   dd_sims[[mc]] <- DAISIE_sim(time = time,
#                               M = M,
#                               pars = init_dd_pars,
#                               replicates = 1,
#                               plot_sims = FALSE,
#                               verbose = FALSE)
# }
#
#
# #Step 7
# di_ml_est_dd <- list()
# dd_ml_est_dd <- list()
# likelihood_ratio_dd <- c()
#
# cat('\nPerforming bootstrap to determine power ...\n')
# for(mc in 1:endmc) {
#   cat('\nAnalyzing simulation:',mc,'\n')
#   di_ml_est_dd[[mc]] <- DAISIE_ML_CS(datalist = dd_sims[[mc]],
#                                      datatype = datatype,
#                                      initparsopt = init_di_pars,
#                                      idparsopt = c(1, 2, 4, 5),
#                                      parsfix = Inf,
#                                      idparsfix = 3,
#                                      idparsnoshift = idparsnoshift,
#                                      idparsmat = idparsmat,
#                                      res = res,
#                                      ddmodel = ddmodel, #should this be set to zero or does it not matter when K is set to inf?
#                                      cond = cond,
#                                      island_ontogeny = island_ontogeny,
#                                      eqmodel = eqmodel,
#                                      x_E = x_E,
#                                      x_I = x_I,
#                                      tol = tol,
#                                      maxiter = maxiter,
#                                      methode = methode,
#                                      optimmethod = optimmethod,
#                                      CS_version = CS_version,
#                                      verbose = verbose,
#                                      tolint = tolint)
#
#   dd_ml_est_dd[[mc]] <- DAISIE_ML_CS(datalist = dd_sims[[mc]],
#                                      datatype = datatype,
#                                      initparsopt = init_dd_pars,
#                                      idparsopt = 1:5,
#                                      parsfix = NULL,
#                                      idparsfix = NULL,
#                                      idparsnoshift = idparsnoshift,
#                                      idparsmat = idparsmat,
#                                      res = res,
#                                      ddmodel = ddmodel,
#                                      cond = cond,
#                                      island_ontogeny = island_ontogeny,
#                                      eqmodel = eqmodel,
#                                      x_E = x_E,
#                                      x_I = x_I,tol = tol,
#                                      maxiter = maxiter,
#                                      methode = methode,
#                                      optimmethod = optimmethod,
#                                      CS_version = CS_version,
#                                      verbose = verbose,
#                                      tolint = tolint)
#   likelihood_ratio[mc] <- dd_ml_est_di$loglik - di_ml_est_di$loglik
# }
# #Step 8
# power <- calc_power(out$LR[(endmc + 2):(2 * endmc + 1)], LRalpha)
#
# funpoweroftest <- function(samples, x) {
#   samplessort <- sort(samples)
#   pup <- which(samplessort > x)
#   power <- length(pup)/(length(samples) + 1)
#   return(power)
# }
