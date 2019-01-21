<<<<<<< HEAD
# # # # # ontogeny test
# # # #
# utils::data(Galapagos_datalist); datalist <- Galapagos_datalist
# Galapagos_datalist
# # # #
# # # # # - pars1[1:4] = Apars
# # # # # - pars1[5] = lac = (initial) cladogenesis rate
# # # # # - pars1[6:7] = extinction rate parameters
# # # # # - pars1[8] = K = maximum number of species possible in the clade
# # # # # - pars1[9] = gam = (initial) immigration rate
# # # # # - pars1[10] = laa = (initial) anagenesis rate
# # # #
# # # # # initparsopt <- c(1000, 0.2, 1, 10, 0.5, 0.05, 100, 0.5, 0.001, 0.2)
# # # # initparsopt <- c(0.5, 0.05, 100, 0.5, 0.001, 0.2)
# # # # idparsfix <- c(1:4)
# # # # parsfix <- c(1000, 0.2, 1, 10)
# # # # idparsopt <- c(5:10)
# # # # ML_out <- DAISIE_ML(
# # # #   datalist = datalist,
# # # #   initparsopt = initparsopt,
# # # #   idparsopt = idparsopt,
# # # #   parsfix = parsfix,
# # # #   idparsfix = idparsfix,
# # # #   island_ontogeny = translate_island_ontogeny("quadratic")
# # # # )
# # # #
# # # #
# # # # # instand_ontogeny = constant test
# # # #
# # # # # initparsopt <- c(1000, 0.2, 1, 10, 0.5, 0.05, 100, 0.5, 0.001, 0.2)
# # # # initparsopt <- c(0.5, 0.05, Inf, 0.001, 0.2)  # mu_min = mu_max for const area
# # # # idparsfix <- c(1:4, 7)
# # # # parsfix <- c(1, 0.2, 1, 10, 0.05)
# # # # idparsopt <- c(5, 6, 8:10)
# # # # ML_out <- DAISIE_ML(
# # # #   datalist = datalist,
# # # #   initparsopt = initparsopt,
# # # #   idparsopt = idparsopt,
# # # #   parsfix = parsfix,
# # # #   idparsfix = idparsfix,
# # # #   island_ontogeny = translate_island_ontogeny("const"),
# # # #   verbose = TRUE
# # # # )
# # # #
# # # # # Old routine
# # # #
=======
# # # # # # # ontogeny test
# # # # # #
# # # # # # utils::data(Galapagos_datalist); datalist <- Galapagos_datalist
# # # # # # Galapagos_datalist
# # # # # #
# # # # # # # - pars1[1:4] = Apars
# # # # # # # - pars1[5] = lac = (initial) cladogenesis rate
# # # # # # # - pars1[6:7] = extinction rate parameters
# # # # # # # - pars1[8] = K = maximum number of species possible in the clade
# # # # # # # - pars1[9] = gam = (initial) immigration rate
# # # # # # # - pars1[10] = laa = (initial) anagenesis rate
# # # # # #
# # # # # # # initparsopt <- c(1000, 0.2, 1, 10, 0.5, 0.05, 100, 0.5, 0.001, 0.2)
# # # # # # initparsopt <- c(0.5, 0.05, 100, 0.5, 0.001, 0.2)
# # # # # # idparsfix <- c(1:4)
# # # # # # parsfix <- c(1000, 0.2, 1, 10)
# # # # # # idparsopt <- c(5:10)
# # # # # # ML_out <- DAISIE_ML(
# # # # # #   datalist = datalist,
# # # # # #   initparsopt = initparsopt,
# # # # # #   idparsopt = idparsopt,
# # # # # #   parsfix = parsfix,
# # # # # #   idparsfix = idparsfix,
# # # # # #   island_ontogeny = translate_island_ontogeny("quadratic")
# # # # # # )
# # # # # #
# # # # # #
# # # # # # # instand_ontogeny = constant test
# # # # # #
# # # # # # # initparsopt <- c(1000, 0.2, 1, 10, 0.5, 0.05, 100, 0.5, 0.001, 0.2)
# initparsopt <- c(0.5, 0.05, Inf, 0.001, 0.2)  # mu_min = mu_max for const area
# idparsfix <- c(1:4, 7)
# parsfix <- c(1, 0.2, 1, 10, 0.05)
# idparsopt <- c(5, 6, 8:10)
# ML_out <- DAISIE_ML(
#   datalist = datalist,
#   initparsopt = initparsopt,
#   idparsopt = idparsopt,
#   parsfix = parsfix,
#   idparsfix = idparsfix,
#   island_ontogeny = translate_island_ontogeny("const"),
#   verbose = TRUE, ddmodel = 11
# )
# # # 
# # # # Old routine
# # # 
>>>>>>> 5bd2e4c06697149bfc4ffca3071a7bd1bb2be089
# initparsopt_cr <- c(0.5, 0.05, Inf, 0.001, 0.2)  # mu_min = mu_max for const area
# idparsfix_cr <- NULL
# parsfix_cr <- NULL
# idparsopt_cr <- c(1:5)
<<<<<<< HEAD
# 
# 
# 
# ML_out_cr <- DAISIE_ML(
#   datalist = datalist,
#   initparsopt = initparsopt_cr,
#   idparsopt = idparsopt_cr,
#   parsfix = parsfix_cr,
#   idparsfix = idparsfix_cr,
#   island_ontogeny = NA,
#   verbose = TRUE,
#   ddmodel = 0
# )
# # # #
# 
# #   island_ontogeny = NA,
# #   verbose = TRUE,
# #   ddmodel = 0
# # )
# # # #
# # # #
# # # #
# # # # # Check specific step (line 33 of log)
# # # #
# # # # pars1_cr_LL = c(0.500000, 1.210526, Inf, 0.001000, 0.200000)
# # # # pars2 = c(40,11,0,0)
# # # # loglik_CS = DAISIE_loglik_all(pars1 = pars1_cr_LL,
# # # #                               pars2 = pars2,
# # # #                               datalist = Galapagos_datalist,
# # # #                               methode = 'ode45')
# # # # pars1_td <- c(max_area = 1,
# # # #               proportional_peak_t = 0.2,
# # # #               peak_sharpness = 1,
# # # #               total_island_age = 15,
# # # #               lac = pars1[1],
# # # #               mu_min = pars1[2],
# # # #               mu_max = pars1[2],
# # # #               K0 = pars1[3],
# # # #               gam = pars1[4],
# # # #               laa = pars1[5])
# # # pars1_td <- order_pars1(pars1_td)
# # # pars2 <- c(pars2,translate_island_ontogeny('const'))
# # # loglik_time <- DAISIE_loglik_all(pars1 = pars1_td,
# # #                                  pars2 = pars2,
# # #                                  datalist = Galapagos_datalist,
# # #                                  methode = "ode45")
# # #
# # #
# # Ks test
# # # ont code
# # test_pars <- c(1, 0.2, 1, 10, 1.142857, 0.329114, 0.329114, 1, 0.251563, 0.578947)
# # pars2 <- c(40,0,0,1);   pars2 <- c(pars2,translate_island_ontogeny('const'))
# # DAISIE_loglik_all(pars1 = test_pars,
# #                   pars2 = pars2,
# #                   datalist = Galapagos_datalist,
# #                   methode = "ode45")
# # 
# # # no ont code
# # no_ont_test_pars <- c(test_pars[5:6], test_pars[8:10])
# # pars2 <- c(40,0,0,1);   pars2 <- c(pars2)
# # DAISIE_loglik_all(pars1 = no_ont_test_pars,
# #                   pars2 = pars2,
# #                   datalist = Galapagos_datalist,
# #                   methode = "ode45")
# 
# 
=======
# 
# 
# ML_out_cr <- DAISIE_ML(
#   datalist = datalist,
#   initparsopt = initparsopt_cr,
#   idparsopt = idparsopt_cr,
#   parsfix = parsfix_cr,
#   idparsfix = idparsfix_cr,
#   island_ontogeny = NA,
#   verbose = TRUE,
#   ddmodel = 11
# )
# # # # # #
# # # # # #
# # # # # #
# # # # # # # Check specific step (line 33 of log)
# # # # # #
# # # # # # pars1_cr_LL = c(0.500000, 1.210526, Inf, 0.001000, 0.200000)
# # # # # # pars2 = c(40,11,0,0)
# # # # # # loglik_CS = DAISIE_loglik_all(pars1 = pars1_cr_LL,
# # # # # #                               pars2 = pars2,
# # # # # #                               datalist = Galapagos_datalist,
# # # # # #                               methode = 'ode45')
# # # # # # pars1_td <- c(max_area = 1,
# # # # # #               proportional_peak_t = 0.2,
# # # # # #               peak_sharpness = 1,
# # # # # #               total_island_age = 15,
# # # # # #               lac = pars1[1],
# # # # # #               mu_min = pars1[2],
# # # # # #               mu_max = pars1[2],
# # # # # #               K0 = pars1[3],
# # # # # #               gam = pars1[4],
# # # # # #               laa = pars1[5])
# # # # # pars1_td <- order_pars1(pars1_td)
# # # # # pars2 <- c(pars2,translate_island_ontogeny('const'))
# # # # # loglik_time <- DAISIE_loglik_all(pars1 = pars1_td,
# # # # #                                  pars2 = pars2,
# # # # #                                  datalist = Galapagos_datalist,
# # # # #                                  methode = "ode45")
# # # # #
# # # # #
# # # Ks test
# # # ont code
# # test_pars <- c(1, 0.2, 1, 10, 1.142857, 0.329114, 0.329114, 1, 0.251563, 0.578947)
# # pars2 <- c(40,11,0,1);   pars2 <- c(pars2,translate_island_ontogeny('const'))
# # DAISIE_loglik_all(pars1 = test_pars,
# #                   pars2 = pars2,
# #                   datalist = Galapagos_datalist,
# #                   methode = "ode45")
# # 
# # # no ont code
# # no_ont_test_pars <- c(test_pars[5:6], test_pars[8:10])
# # pars2 <- c(40,11,0,1);   pars2 <- c(pars2)
# # DAISIE_loglik_all(pars1 = no_ont_test_pars,
# #                   pars2 = pars2,
# #                   datalist = Galapagos_datalist,
# #                   methode = "ode45")
# # 
# # 
>>>>>>> 5bd2e4c06697149bfc4ffca3071a7bd1bb2be089
