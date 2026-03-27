# The expected number of endemics and non-endemics under the DAISIE model

This function calculates the expected number of endemics, non-endemics
and the sum of these for a given set of parameter values, a given
mainland species pool size and a given time, where there can be
diversity-dependence

## Usage

``` r
DAISIE_ExpEIN2(
  tvec,
  pars,
  M,
  initEI = NULL,
  res = 1000,
  ddmodel = 11,
  methode = "ode45",
  reltolint = 1e-16,
  abstolint = 1e-16
)
```

## Arguments

- tvec:

  The times at which the probabilities need to be computed.

- pars:

  A numeric vector containing the model parameters:

  - `pars[1]`: lambda^c (cladogenesis rate)

  - `pars[2]`: mu (extinction rate)

  - `pars[3]`: K (carrying capacity), set K=Inf for diversity
    independence.

  - `pars[4]`: gamma (immigration rate)

  - `pars[5]`: lambda^a (anagenesis rate)

  - `pars[6]`: lambda^c (cladogenesis rate) for either type 2 species or
    rate set 2 in rate shift model

  - `pars[7]`: mu (extinction rate) for either type 2 species or rate
    set 2 in rate shift model

  - `pars[8]`: K (carrying capacity) for either type 2 species or rate
    set 2 in rate shift model, set K=Inf for diversity independence.

  - `pars[9]`: gamma (immigration rate) for either type 2 species or
    rate set 2 in rate shift model

  - `pars[10]`: lambda^a (anagenesis rate) for either type 2 species or
    rate set 2 in rate shift model

  Elements 6:10 are required only when type 2 species are included or in
  the rate shift model. For
  [`DAISIE_sim_relaxed_rate()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim_relaxed_rate.md)
  `pars[6]` is the standard deviation of the gamma distribution for the
  relaxed parameter and the parameter chosen by the `relaxed_par`
  argument is the mean of the gamma distribution for the relaxed
  parameter.

- M:

  Numeric defining the size of mainland pool, i.e. the number of species
  that can potentially colonize the island.

- initEI:

  The initial values for the number of endemics and non-endemics. In
  [`DAISIE_probdist()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_probdist.md)
  or
  [`DAISIE_margprobdist()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_margprobdist.md)
  either this or initprobs must be NULL. In
  [`DAISIE_numcol()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_numcol.md)
  when it is NULL, it is assumed that the island is empty.

- res:

  Sets the maximum number of species for which a probability must be
  computed, must be larger than the size of the largest clade.

- ddmodel:

  Sets the model of diversity-dependence:  
    

  - ddmodel = 0 : no diversity dependence

  - ddmodel = 1 : linear dependence in speciation rate

  - ddmodel = 11: linear dependence in speciation rate and in
    immigration rate

  - ddmodel = 2 : exponential dependence in speciation rate

  - ddmodel = 21: exponential dependence in speciation rate and in
    immigration rate

- methode:

  Method of the ODE-solver. Supported Boost `ODEINT` solvers (steppers)
  are: `"odeint::runge_kutta_cash_karp54"`
  `"odeint::runge_kutta_fehlberg78"` `"odeint::runge_kutta_dopri5"`
  `"odeint::bulirsch_stoer"` without `odeint::`-prefix,
  [`ode`](https://rdrr.io/pkg/deSolve/man/ode.html) method is assumed.
  The default method overall is `"lsodes"` for
  [`DAISIE_ML_CS()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_ML.md)
  and `"ode45"` from [`ode()`](https://rdrr.io/pkg/deSolve/man/ode.html)
  for
  [`DAISIE_ML_IW()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_ML_IW.md).

- reltolint:

  Numeric relative tolerance of the integration

- abstolint:

  Numeric absolute tolerance of the integration

## Value

- tot_expEIN:

  The output is a list with three elements:  
    
  `ExpE` The number of endemic species at the times in tvec  
  `ExpI` The number of non-endemic species at the times in tvec  
  `ExpN` The sum of the number of endemics and non-endemics at the times
  in tvec

## Author

Rampal S. Etienne

## Examples

``` r
DAISIE_ExpEIN2(tvec = c(0.000001,0.5,0.75,1),
                         pars = c(0.3,0.1,10,1,0.1),
                         M = 1000,
                         initEI = rbind(c(1,0),c(2,0),c(0,1)))
#> $ExpE
#> [1]   3.000001  68.398821 136.540451 219.906946
#> 
#> $ExpI
#> [1]   1.000998 353.637331 452.384821 519.722809
#> 
#> $ExpN
#> [1]   4.000999 422.036152 588.925273 739.629755
#> 
```
