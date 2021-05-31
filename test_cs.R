library(DAISIE)

pars1 <- c(2.000, 2.700, 20.000, 0.009, 1.010)
pars2 <- c(1.0e+02, 1.1e+01, 0.0e+00, 0.0e+00, NA, 0.0e+00, 1.0e-04,
           1.0e-05, 1.0e-07, 3.0e+03, 9.5e-01, 9.8e-01)
brts <- c(4.0000, 3.0282, 1.3227, 0.8223, 0.4286, 0.3462, 0.2450, 0.0808,
          0.0527, 0.0327, 0.0221, 0.1180, 0.0756, 0.0525, 0.0322, 0.0118)
stac <- 2
missnumspec <- 0
loglik1 <- DAISIE_loglik_CS_choice(pars1 = pars1,
                                  pars2 = pars2,
                                  brts = brts,
                                  stac = stac,
                                  missnumspec = missnumspec)
loglik2 <- DAISIE_loglik_CS_choice(pars1 = pars1,
                                  pars2 = pars2,
                                  brts = brts,
                                  stac = stac,
                                  missnumspec = missnumspec,
                                  methode = "odeint::runge_kutta_cash_karp54")
print(c(loglik1, loglik2, -17.6550433826))



# # deSolve: warnings, odeint: too many steps
# pars1 <- c(2.000, 2.700, 20.000, 0.009, 1.010)
# pars2 <- c(1.0e+02, 1.1e+01, 0.0e+00, 0.0e+00, NA, 0.0e+00, 1.0e-04,
#            1.0e-05, 1.0e-07, 3.0e+03, 9.5e-01, 9.8e-01)
# brts <- c(4.0000, 3.0282, 1.3227, 0.8223, 0.4286, 0.3462, 0.2450, 0.0808,
#           0.0527, 0.0327, 0.0221, 0.1180, 0.0756, 0.0525, 0.0322, 0.0118)
# stac <- 2
# missnumspec <- 0
# CS_version <- list(model = 2,
#                    relaxed_par = "cladogenesis",
#                    sd = 1)
# loglik1 <- DAISIE_loglik_CS_choice(pars1 = pars1,
#                                   pars2 = pars2,
#                                   brts = brts,
#                                   stac = stac,
#                                   missnumspec = missnumspec,
#                                   CS_version = CS_version)
# loglik2 <- DAISIE_loglik_CS_choice(pars1 = pars1,
#                                    pars2 = pars2,
#                                    brts = brts,
#                                    stac = stac,
#                                    missnumspec = missnumspec,
#                                    CS_version = CS_version,
#                                    methode = "odeint::runge_kutta_fehlberg78")
# print(c(loglik1, loglik2, -9.55117524011))


pars1 <- c(2.000, 2.700, 20.000, 0.009, 1.010)
pars2 <- c(1.0e+02, 1.1e+01, 0.0e+00, 0.0e+00, NA, 0.0e+00, 1.0e-04,
           1.0e-05, 1.0e-07, 3.0e+03, 9.5e-01, 9.8e-01)
brts <- c(4.0000, 3.0282, 1.3227, 0.8223, 0.4286, 0.3462, 0.2450, 0.0808,
          0.0527, 0.0327, 0.0221, 0.1180, 0.0756, 0.0525, 0.0322, 0.0118)
stac <- 2
missnumspec <- 0
CS_version <- 0
loglik1 <- DAISIE_loglik_CS_choice(pars1 = pars1,
                                   pars2 = pars2,
                                   brts = brts,
                                   stac = stac,
                                   missnumspec = missnumspec,
                                   CS_version = CS_version)
loglik2 <- DAISIE_loglik_CS_choice(pars1 = pars1,
                                   pars2 = pars2,
                                   brts = brts,
                                   stac = stac,
                                   missnumspec = missnumspec,
                                   CS_version = CS_version,
                                   methode = "odeint::runge_kutta_fehlberg78")
print(c(loglik1, loglik2))



# # deSolve: roundoff error, odeint: too many steps
# utils::data(Galapagos_datalist)
# loglik1 <- DAISIE::DAISIE_loglik_all(pars1 = c(2.55068735, 2.68345455, 10.00000000, 0.00933207, 1.01007312),
#                                      pars2 = c(100, 0, 0, 0, NA),
#                                      datalist = Galapagos_datalist,
#                                      CS_version = list(model = 2,
#                                                        relaxed_par = "cladogenesis",
#                                                        sd = 1),
#                                      abstolint = 1e-18,
#                                      reltolint = 1e-20)
# loglik2 <- DAISIE::DAISIE_loglik_all(pars1 = c(2.55068735, 2.68345455, 10.00000000, 0.00933207, 1.01007312),
#                                      pars2 = c(100, 0, 0, 0, NA),
#                                      datalist = Galapagos_datalist,
#                                      methode = "odeint::runge_kutta_fehlberg78",
#                                      CS_version = list(model = 2,
#                                                        relaxed_par = "cladogenesis",
#                                                        sd = 1),
#                                      abstolint = 1e-18,
#                                      reltolint = 1e-20)
# print(c(loglik1, loglik2, --77.5108137039949))

