test_that("various solvers give the same result for max age unequal to island age", {
  data("NewZealand_birds_datalist")
  datalist <- NewZealand_birds_datalist
  i <- 35
  pars1 = c(0.1255591, 0.12, Inf, 0.4697986, 0.326698)
  pars2 <- c(25, 11,0,2)
  loglik_DAISIE_DE <- DAISIE_DE_logpES_max_age_coltime(datalist, i, pars1 = c(pars1[1:2],pars1[2],pars1[4:5]), methode = 'lsodes', atol = 1E-16, rtol = 1E-16)
  loglik_DAISIE_R <- DAISIE:::DAISIE_loglik_CS_choice(pars1 = pars1,
                                         pars2 = pars2,
                                         datalist = datalist,
                                         brts = datalist[[i]]$branching_times,
                                         stac = datalist[[i]]$stac,
                                         missnumspec = datalist[[i]]$missing_species,
                                         methode = "deSolve_R::lsodes")
  loglik_DAISIE_FORTRAN <- DAISIE:::DAISIE_loglik_CS_choice(pars1 = pars1,
                                         pars2 = pars2,
                                         datalist = datalist,
                                         brts = datalist[[i]]$branching_times,
                                         stac = datalist[[i]]$stac,
                                         missnumspec = datalist[[i]]$missing_species,
                                         methode = "lsodes")
  testthat::expect_equal(loglik_DAISIE_DE,loglik_DAISIE_R, tol = 1E-3)
  testthat::expect_equal(loglik_DAISIE_FORTRAN,loglik_DAISIE_R)
  pars1 = c(0.1255591, 0.12, 50, 0.4697986, 0.326698)
  pars2 <- c(25, 11,0,2)
  loglik_DAISIE_R <- DAISIE:::DAISIE_loglik_CS_choice(pars1 = pars1,
                                                      pars2 = pars2,
                                                      datalist = datalist,
                                                      brts = datalist[[i]]$branching_times,
                                                      stac = datalist[[i]]$stac,
                                                      missnumspec = datalist[[i]]$missing_species,
                                                      methode = "deSolve_R::lsodes")
  loglik_DAISIE_FORTRAN <- DAISIE:::DAISIE_loglik_CS_choice(pars1 = pars1,
                                                            pars2 = pars2,
                                                            datalist = datalist,
                                                            brts = datalist[[i]]$branching_times,
                                                            stac = datalist[[i]]$stac,
                                                            missnumspec = datalist[[i]]$missing_species,
                                                            methode = "lsodes")
  testthat::expect_equal(loglik_DAISIE_FORTRAN,loglik_DAISIE_R)
})
