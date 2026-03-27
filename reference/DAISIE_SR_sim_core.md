# DEPRECATED - Algorithm component of DAISIE_SR_sim.

\#' @details This function's use has been deprecated in favour of
[`DAISIE_sim_cr_shift()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim_cr_shift.md).
Please use that function instead.

## Usage

``` r
DAISIE_SR_sim_core(time, mainland_n, pars)
```

## Arguments

- time:

  Length of the simulation in time units. For example, if an island is
  know to be 4 million years old, setting time = 4 will simulate entire
  life span of the island; setting time = 2 will stop the simulation at
  the mid-life of the island.

- mainland_n:

  A numeric stating the number of mainland species, that is the number
  of species that can potentially colonize the island. If using a
  clade-specific diversity dependence, this value is set to 1. If using
  an island-wide diversity dependence, this value is set to the number
  of mainland species.

- pars:

  Contains the model parameters:  
    
  `pars[1]` corresponds to lambda^c (cladogenesis rate) before the
  shift  
  `pars[2]` corresponds to mu (extinction rate) before the shift  
  `pars[3]` corresponds to K (clade-level carrying capacity) before the
  shift. Set K=Inf for non-diversity dependence.  
  `pars[4]` corresponds to gamma (immigration rate) before the shift  
  `pars[5]` corresponds to lambda^a (anagenesis rate) before the shift  
  `pars[6]` corresponds to lambda^c (cladogenesis rate) after the
  shift  
  `pars[7]` corresponds to mu (extinction rate) after the shift  
  `pars[8]` corresponds to K (clade-level carrying capacity) after the
  shift. Set K=Inf for non-diversity dependence.  
  `pars[9]` corresponds to gamma (immigration rate) after the shift  
  `pars[10]` corresponds to lambda^a (anagenesis rate) after the shift  
  `pars[11]` corresponds to the time of shift. This is defined as time
  before the end of the simulation. For example, setting time = 4 and
  pars\[11\] = 1.5 will simulate with pars\[1:5\] from 4 to 1.5 and with
  pars\[6:10\] from 1.5 to 0.

## Value

List with DAISIE simulation.

## References

Hauffe, T., D. Delicado, R.S. Etienne and L. Valente (2020). Lake
expansion increases equilibrium diversity via the target effect of
island biogeography

## Author

Luis Valente, Albert Phillimore, and Torsten Hauffe
