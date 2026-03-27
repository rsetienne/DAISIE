# Internal function of the DAISIE simulation

Internal function of the DAISIE simulation

## Usage

``` r
DAISIE_sim_core_cr(
  time,
  mainland_n,
  pars,
  nonoceanic_pars,
  hyper_pars,
  area_pars,
  use_rcpp = FALSE
)
```

## Arguments

- time:

  Numeric defining the length of the simulation in time units. For
  example, if an island is known to be 4 million years old, setting time
  = 4 will simulate the entire life span of the island; setting time = 2
  will stop the simulation at the mid-life of the island.

- mainland_n:

  A numeric stating the number of mainland species, that is the number
  of species that can potentially colonize the island. If using a
  clade-specific diversity dependence, this value is set to 1. If using
  an island-wide diversity dependence, this value is set to the number
  of mainland species.

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

- nonoceanic_pars:

  A vector of length two with:

  - \[1\]: the probability of sampling a species from the mainland

  - \[2\]: the probability of the species sampled from the mainland
    being nonendemic

- hyper_pars:

  A named list of numeric hyperparameters for the rate calculations as
  returned by
  [`create_hyper_pars()`](https://rsetienne.github.io/DAISIE/reference/create_hyper_pars.md):

  - \[1\]: is d the scaling parameter for exponent for calculating
    cladogenesis rate

  - \[2\]: is x the exponent for calculating extinction rate

- area_pars:

  A named list containing area and sea level parameters as created by
  [`create_area_pars()`](https://rsetienne.github.io/DAISIE/reference/create_area_pars.md):

  - \[1\]: maximum area

  - \[2\]: current area

  - \[3\]: value from 0 to 1 indicating where in the island's history
    the peak area is achieved

  - \[4\]: total island age

  - \[5\]: amplitude of area fluctuation from sea level

  - \[6\]: frequency of sine wave of area change from sea level

  - \[7\]: angle of the slope of the island

- use_rcpp:

  If TRUE, use Rcpp implementation of DAISIE simulation core. Default is
  FALSE.
