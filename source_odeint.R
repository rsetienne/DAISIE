# source_odeint.R
#
# partial odeint test suite from tests/testthat/test_odeint.R
# Run with:
#
# R -d "valgrind --tool=memcheck --leak-check=full --error-exitcode=1" --vanilla -e "source('source_odeint.R')" &> source_odeint.log


library(DAISIE)

utils::data(Galapagos_datalist_2types)
pars1 <- c(
0.195442017,
0.087959583,
Inf,
0.002247364,
0.873605049,
3755.202241,
8.909285094,
14.99999923,
0.002247364,
0.873605049,
0.163
)
pars2 <- c(40, 11, 0, 1)
methode <- 'deSolve_R::lsodes'
loglik_lsodes_R <- DAISIE_loglik_all(pars1, pars2, Galapagos_datalist_2types, methode = methode)
methode <- 'deSolve_R::lsoda'
loglik_lsoda_R <- DAISIE_loglik_all(pars1, pars2, Galapagos_datalist_2types, methode = methode)
methode <- 'lsodes'
loglik_lsodes <- DAISIE_loglik_all(pars1, pars2, Galapagos_datalist_2types, methode = methode)
methode <- 'lsoda'
loglik_lsoda <- DAISIE_loglik_all(pars1, pars2, Galapagos_datalist_2types, methode = methode)
methode <- 'odeint::runge_kutta_cash_karp54'
loglik_rkck54 <- DAISIE_loglik_all(pars1, pars2, Galapagos_datalist_2types, methode = methode)
methode <- 'odeint::runge_kutta_fehlberg78'
loglik_rkf78 <- DAISIE_loglik_all(pars1, pars2, Galapagos_datalist_2types, methode = methode)
methode <- 'odeint::runge_kutta_dopri5'
loglik_rkd5 <- DAISIE_loglik_all(pars1, pars2, Galapagos_datalist_2types, methode = methode)
methode <- 'odeint::bulirsch_stoer'
loglik_bs <- DAISIE_loglik_all(pars1, pars2, Galapagos_datalist_2types, methode = methode)
methode <- 'odeint::rosenbrock4'
loglik_rb <- DAISIE_loglik_all(pars1, pars2, Galapagos_datalist_2types, methode = methode)
methode <- 'odeint::adams_bashforth_moulton_1'
DAISIE_CS_max_steps(100000000)
DAISIE_abm_factor(0.000001)
loglik_abm <- DAISIE_loglik_all(pars1, pars2, Galapagos_datalist_2types, methode = methode)
