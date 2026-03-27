# Simulate (non-)oceanic islands with given parameters under a rate-shift regime

This function simulates islands with given cladogenesis, extinction,
Kprime, immigration and anagenesis parameters, all of which modelled as
time-constant parameters, which can be switched to a different
diversification regime (i.e., different set of parameters) at one or
more set times before the present. Further, it allows for the simulation
of non-oceanic islands, generating islands for which the starting
condition includes potential endemic and non-endemic species.

## Usage

``` r
DAISIE_sim_cr_shift(
  time,
  M,
  pars,
  replicates,
  shift_times,
  divdepmodel = "CS",
  nonoceanic_pars = c(0, 0),
  num_guilds = NULL,
  sample_freq = 25,
  plot_sims = TRUE,
  hyper_pars = create_hyper_pars(d = 0, x = 0),
  area_pars = DAISIE::create_area_pars(max_area = 1, current_area = 1,
    proportional_peak_t = 0, total_island_age = 0, sea_level_amplitude = 0,
    sea_level_frequency = 0, island_gradient_angle = 0),
  cond = 0,
  verbose = TRUE,
  ...
)
```

## Arguments

- time:

  Numeric defining the length of the simulation in time units. For
  example, if an island is known to be 4 million years old, setting time
  = 4 will simulate the entire life span of the island; setting time = 2
  will stop the simulation at the mid-life of the island.

- M:

  Numeric defining the size of mainland pool, i.e. the number of species
  that can potentially colonize the island.

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

- replicates:

  Integer specifying number of island replicates to be simulated.

- shift_times:

  a numeric vector specifying when the rate shifts occur before the
  present.

- divdepmodel:

  Option divdepmodel = 'CS' runs a model with clade-specific carrying
  capacity, where diversity-dependence operates only within single
  clades, i.e. only among species originating from the same mainland
  colonist. Option divdepmodel = 'IW' runs a model with island-wide
  carrying capacity, where diversity-dependence operates within and
  among clades. Option divdepmodel = 'GW' runs a model with
  diversity-dependence operates within a guild.

- nonoceanic_pars:

  A vector of length two with:

  - \[1\]: the probability of sampling a species from the mainland

  - \[2\]: the probability of the species sampled from the mainland
    being nonendemic

- num_guilds:

  The number of guilds on the mainland. The number of mainland species
  is divided by the number of guilds when `divdepmodel = "GW"`

- sample_freq:

  Numeric specifing the number of units times should be divided by for
  plotting purposes. Larger values will lead to plots with higher
  resolution, but will also run slower.

- plot_sims:

  `Default = TRUE` plots species-through-time (STT) plots. It detects
  how many types of species are present. If only one type of species is
  present, STT is plotted for all species. If two types are present,
  three plots are produced: STT for all, STT for type 1 and STT for type
  2.

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

- cond:

  cond = 0 : conditioning on island age  
  cond = 1 : conditioning on island age and non-extinction of the island
  biota  
  . cond \> 1 : conditioning on island age and having at least cond
  colonizations on the island. This last option is not yet available for
  the IW model  

- verbose:

  A numeric vector of length 1, which in simulations and
  \`DAISIEdataprep()\` can be \`1\` or \`0\`, where \`1\` gives
  intermediate output should be printed. For ML functions a numeric
  determining if intermediate output should be printed. The default:
  \`0\` does not print, \`1\` prints the initial likelihood and the
  settings that were selected (which parameters are to be optimised,
  fixed or shifted), \`2\` prints the same as \`1 and also the
  intermediate output of the parameters and loglikelihood, while \`3\`
  the same as \`2\` and prints intermediate progress during likelihood
  computation.

- ...:

  Any arguments to pass on to plotting functions.

## Value

A list. The highest level of the least corresponds to each individual
replciate. The first element of each replicate is composed of island
information containing:

- `$island_age`: A numeric with the island age.

- `$not_present`: A numeric with the number of mainland lineages that
  are not present on the island.

- `$stt_all`: STT table for all species on the island (nI - number of
  non-endemic species; nA - number of anagenetic species, nC - number of
  cladogenetic species, present - number of independent colonisations
  present)

- `$brts_table`: Only for simulations under `"IW"`. Table containing
  information on order of events in the data, for use in maximum
  likelihood optimization.).

The subsequent elements of the list pertaining to each replcate contain
information on a single colonist lineage on the island and have 4
components:

- `$branching_times`: island age and stem age of the population/species
  in the case of Non-endemic, Non-endemic_MaxAge and Endemic anagenetic
  species. For cladogenetic species these should be island age and
  branching times of the radiation including the stem age of the
  radiation.

- `$stac`: An integer ranging from 1 to 4 indicating the status of the
  colonist:

  1.  Non_endemic_MaxAge

  2.  Endemic

  3.  Endemic&Non_Endemic

  4.  Non_endemic_MaxAge

- `$missing_species`: number of island species that were not sampled for
  particular clade (only applicable for endemic clades)

- `$type_1or2`: whether the colonist belongs to type 1 or type 2

## References

Valente, L.M., A.B. Phillimore and R.S. Etienne (2015). Equilibrium and
non-equilibrium dynamics simultaneously operate in the Galapagos
islands. Ecology Letters 18: 844-852.

Hauffe, T., D. Delicado, R.S. Etienne and L. Valente (2020). Lake
expansion elevates equilibrium diversity via increasing colonization.

## See also

[`DAISIE_plot_sims()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_plot_sims.md)
for plotting STT of simulation outputs.

Other simulation models:
[`DAISIE_sim()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim.md),
[`DAISIE_sim_relaxed_rate()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim_relaxed_rate.md),
[`DAISIE_sim_time_dep()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim_time_dep.md)

## Author

Luis Valente, Albert Phillimore, Torsten Hauffe
