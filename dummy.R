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

methode <- 'odeint::runge_kutta_cash_karp54'
loglik_rkck54 <- DAISIE_loglik_all(pars1, pars2, Galapagos_datalist_2types, methode = methode)

methode <- 'odeint::bulirsch_stoer'
loglik_bs <- DAISIE_loglik_all(pars1, pars2, Galapagos_datalist_2types, methode = methode)
