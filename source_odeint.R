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


==2373== Memcheck, a memory error detector
==2373== Copyright (C) 2002-2022, and GNU GPL'd, by Julian Seward et al.
==2373== Using Valgrind-3.19.0 and LibVEX; rerun with -h for copyright info
==2373== Command: /home/hanno/opt/Rval/lib/R/bin/exec/R --vanilla -e source('source_odeint.R')
==2373== 

R Under development (unstable) (2023-03-07 r83950) -- "Unsuffered Consequences"
Copyright (C) 2023 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> source('source_odeint.R')
Status of colonist: 0, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -0.007835
Status of colonist: 0, Parameters: 3755.202241, 8.909285, 14.999999, 0.002247, 0.873605, Loglikelihood: -0.008967
Status of colonist: 1, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.259846
Status of colonist: 4, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.497364
Status of colonist: 2, Parameters: 3755.202241, 8.909285, 14.999999, 0.002247, 0.873605, Loglikelihood: 1.545883
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -13.576394
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.966530
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -8.752420
Status of colonist: 1, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.259846
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.981045
Parameters: 
lambda^c, mu, K, gamma, lambda^a
0.195442, 0.087960, Inf, 0.002247, 0.873605
Loglikelihood: -61.702819

Status of colonist: 0, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -0.007835
Status of colonist: 0, Parameters: 3755.202241, 8.909285, 14.999999, 0.002247, 0.873605, Loglikelihood: -0.008967
Status of colonist: 1, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.259846
Status of colonist: 4, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.497364
Status of colonist: 2, Parameters: 3755.202241, 8.909285, 14.999999, 0.002247, 0.873605, Loglikelihood: 1.545883
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -13.576394
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.966530
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -8.752420
Status of colonist: 1, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.259846
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.981045
Parameters: 
lambda^c, mu, K, gamma, lambda^a
0.195442, 0.087960, Inf, 0.002247, 0.873605
Loglikelihood: -61.702819

Status of colonist: 0, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -0.007835
Status of colonist: 0, Parameters: 3755.202241, 8.909285, 14.999999, 0.002247, 0.873605, Loglikelihood: -0.008967
Status of colonist: 1, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.259846
Status of colonist: 4, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.497364
Status of colonist: 2, Parameters: 3755.202241, 8.909285, 14.999999, 0.002247, 0.873605, Loglikelihood: 1.545883
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -13.576394
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.966530
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -8.752420
Status of colonist: 1, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.259846
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.981045
Parameters: 
lambda^c, mu, K, gamma, lambda^a
0.195442, 0.087960, Inf, 0.002247, 0.873605
Loglikelihood: -61.702819

Status of colonist: 0, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -0.007835
Status of colonist: 0, Parameters: 3755.202241, 8.909285, 14.999999, 0.002247, 0.873605, Loglikelihood: -0.008967
Status of colonist: 1, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.259846
Status of colonist: 4, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.497364
Status of colonist: 2, Parameters: 3755.202241, 8.909285, 14.999999, 0.002247, 0.873605, Loglikelihood: 1.545883
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -13.576394
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.966530
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -8.752420
Status of colonist: 1, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.259846
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.981045
Parameters: 
lambda^c, mu, K, gamma, lambda^a
0.195442, 0.087960, Inf, 0.002247, 0.873605
Loglikelihood: -61.702819

Status of colonist: 0, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -0.007835
Status of colonist: 0, Parameters: 3755.202241, 8.909285, 14.999999, 0.002247, 0.873605, Loglikelihood: -0.008967
Status of colonist: 1, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.259846
Status of colonist: 4, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.497364
Status of colonist: 2, Parameters: 3755.202241, 8.909285, 14.999999, 0.002247, 0.873605, Loglikelihood: 1.545883
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -13.576394
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.966530
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -8.752420
Status of colonist: 1, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.259846
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.981045
Parameters: 
lambda^c, mu, K, gamma, lambda^a
0.195442, 0.087960, Inf, 0.002247, 0.873605
Loglikelihood: -61.702819

Status of colonist: 0, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -0.007835
Status of colonist: 0, Parameters: 3755.202241, 8.909285, 14.999999, 0.002247, 0.873605, Loglikelihood: -0.008967
Status of colonist: 1, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.259846
Status of colonist: 4, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.497364
Status of colonist: 2, Parameters: 3755.202241, 8.909285, 14.999999, 0.002247, 0.873605, Loglikelihood: 1.545883
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -13.576394
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.966530
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -8.752420
Status of colonist: 1, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.259846
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.981045
Parameters: 
lambda^c, mu, K, gamma, lambda^a
0.195442, 0.087960, Inf, 0.002247, 0.873605
Loglikelihood: -61.702819

Status of colonist: 0, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -0.007835
Status of colonist: 0, Parameters: 3755.202241, 8.909285, 14.999999, 0.002247, 0.873605, Loglikelihood: -0.008967
Status of colonist: 1, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.259846
Status of colonist: 4, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.497364
Status of colonist: 2, Parameters: 3755.202241, 8.909285, 14.999999, 0.002247, 0.873605, Loglikelihood: 1.545883
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -13.576394
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.966530
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -8.752420
Status of colonist: 1, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.259846
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.981045
Parameters: 
lambda^c, mu, K, gamma, lambda^a
0.195442, 0.087960, Inf, 0.002247, 0.873605
Loglikelihood: -61.702819

Status of colonist: 0, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -0.007835
Status of colonist: 0, Parameters: 3755.202241, 8.909285, 14.999999, 0.002247, 0.873605, Loglikelihood: -0.008967
Status of colonist: 1, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.259846
Status of colonist: 4, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.497364
Status of colonist: 2, Parameters: 3755.202241, 8.909285, 14.999999, 0.002247, 0.873605, Loglikelihood: 1.545883
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -13.576394
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.966530
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -8.752420
Status of colonist: 1, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.259846
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.981045
Parameters: 
lambda^c, mu, K, gamma, lambda^a
0.195442, 0.087960, Inf, 0.002247, 0.873605
Loglikelihood: -61.702819

Status of colonist: 0, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -0.007835
Status of colonist: 0, Parameters: 3755.202241, 8.909285, 14.999999, 0.002247, 0.873605, Loglikelihood: -0.008967
Status of colonist: 1, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.259846
Status of colonist: 4, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.497364
Status of colonist: 2, Parameters: 3755.202241, 8.909285, 14.999999, 0.002247, 0.873605, Loglikelihood: 1.545883
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -13.576394
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.966530
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -8.752420
Status of colonist: 1, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.259846
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.981045
Parameters: 
lambda^c, mu, K, gamma, lambda^a
0.195442, 0.087960, Inf, 0.002247, 0.873605
Loglikelihood: -61.702819

Status of colonist: 0, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -0.007835
Status of colonist: 0, Parameters: 3755.202241, 8.909285, 14.999999, 0.002247, 0.873605, Loglikelihood: -0.008967
Status of colonist: 1, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.259847
Status of colonist: 4, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.497364
Status of colonist: 2, Parameters: 3755.202241, 8.909285, 14.999999, 0.002247, 0.873605, Loglikelihood: 1.545893
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -13.576392
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.966530
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -8.752420
Status of colonist: 1, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.259847
Status of colonist: 2, Parameters: 0.195442, 0.087960, Inf, 0.002247, 0.873605, Loglikelihood: -6.981045
Parameters: 
lambda^c, mu, K, gamma, lambda^a
0.195442, 0.087960, Inf, 0.002247, 0.873605
Loglikelihood: -61.702808

> 
> 
==2373== 
==2373== HEAP SUMMARY:
==2373==     in use at exit: 67,395,855 bytes in 14,084 blocks
==2373==   total heap usage: 75,008,366 allocs, 74,994,282 frees, 40,220,411,885 bytes allocated
==2373== 
==2373== LEAK SUMMARY:
==2373==    definitely lost: 0 bytes in 0 blocks
==2373==    indirectly lost: 0 bytes in 0 blocks
==2373==      possibly lost: 0 bytes in 0 blocks
==2373==    still reachable: 67,395,855 bytes in 14,084 blocks
==2373==                       of which reachable via heuristic:
==2373==                         newarray           : 4,264 bytes in 1 blocks
==2373==         suppressed: 0 bytes in 0 blocks
==2373== Reachable blocks (those to which a pointer was found) are not shown.
==2373== To see them, rerun with: --leak-check=full --show-leak-kinds=all
==2373== 
==2373== For lists of detected and suppressed errors, rerun with: -s
==2373== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 0 from 0)
