# Simulate (non-)oceanic islands with given parameters under time-constant rates

This function simulates islands with given cladogenesis, extinction,
Kprime, immigration and anagenesis parameters, all of which modelled as
time-constant parameters. If a single parameter set is provided (5
parameters) it simulates islands where all species have the same
macro-evolutionary process. If two paramater sets (10 parameters) are
provided, it simulates islands where two different macro-evolutionary
processes operate, one applying to type 1 species and other to type 2
species. Further, it allows for the simulation of non-oceanic islands,
generating islands for which the starting condition includes potential
endemic and non-endemic species.

## Usage

``` r
DAISIE_sim_cr(
  time,
  M,
  pars,
  replicates,
  divdepmodel = "CS",
  nonoceanic_pars = c(0, 0),
  num_guilds = NULL,
  prop_type2_pool = NA,
  replicates_apply_type2 = TRUE,
  sample_freq = 25,
  plot_sims = TRUE,
  hyper_pars = create_hyper_pars(d = 0, x = 0),
  area_pars = create_area_pars(max_area = 1, current_area = 1, proportional_peak_t = 0,
    total_island_age = 0, sea_level_amplitude = 0, sea_level_frequency = 0,
    island_gradient_angle = 0),
  cond = 0,
  verbose = TRUE,
  files_to_write = FALSE,
  use_rcpp = FALSE,
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

- prop_type2_pool:

  Fraction of mainland species that belongs to the second subset of
  species (type 2). Applies only when two types of species are simulated
  (length(pars) = 10). For
  [`DAISIE_dataprep()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_dataprep.md)
  applies only if number_clade_types = 2. In
  [`DAISIE_dataprep()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_dataprep.md)
  the default `"proportional"` sets the fraction to be proportional to
  the number of clades of distinct macroevolutionary process that have
  colonised the island.

- replicates_apply_type2:

  Applies only when two types of species are being simulated. Default
  replicates_apply_type2 = TRUE runs simulations until the number of
  islands where a type 2 species has colonised is equal to the specified
  number of replicates. This is recommended if `prop_type2_pool` is
  small or if the rate of immigration of type two species (`pars[9]`) is
  low, meaning that more replicates are needed to achieved an adequate
  sample size of islands with type 2 species. Setting
  `replicates_apply_type2 = FALSE` simulates islands up to the specified
  number of replicates regardless of whether type 2 species have
  colonised or not.

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

- files_to_write:

  number of files to write simulations to file

- use_rcpp:

  If TRUE, use Rcpp implementation of DAISIE simulation core. Default is
  FALSE.

- ...:

  Any arguments to pass on to plotting functions.

## Value

A list. The highest level of the least corresponds to each individual
replciate. The first element of each replicate is composed of island
information containing:

- `$island_age`: A numeric with the island age.

- `$not_present`: the number of mainland lineages that are not present
  on the island. It is only present if only 1 type of species is
  simulated. Becomes `$not_present_type1`: the number of mainland
  lineages of type 1 that are not present on the island and
  `$not_present_type2`: the number of mainland lineages of type 2 that
  are not present on the island, if two types are simulated.

- `$stt_all`: STT table for all species on the island (nI - number of
  non-endemic species; nA - number of anagenetic species, nC - number of
  cladogenetic species, present - number of independent colonisations
  present)

- `$stt_stt_type1`: STT table for type 1 species on the island - only if
  2 types of species were simulated (nI - number of non-endemic species;
  nA - number of anagenetic species, nC - number of cladogenetic
  species, present - number of independent colonisations present).

- `$stt_stt_type2`: STT table for type 2 species on the island - only if
  2 types of species were simulated (nI - number of non-endemic species;
  nA - number of anagenetic species, nC - number of cladogenetic
  species, present - number of independent colonisations present ).

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

## See also

[`DAISIE_plot_sims()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_plot_sims.md)
for plotting STT of simulation outputs.

Other simulation models:
[`DAISIE_sim_cr_shift()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim_cr_shift.md),
[`DAISIE_sim_relaxed_rate()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim_relaxed_rate.md),
[`DAISIE_sim_time_dep()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim_time_dep.md)

## Author

Luis Valente, Albert Phillimore, Joshua Lambert, Shu Xie, Pedro Neves,
Richèl J. C. Bilderbeek, Rampal Etienne

## Examples

``` r
## Simulate 2 islands for 1 million years, where all species have equal
## rates. Pool size 100.

clado_rate <- 0.5
ext_rate <- 0.2
carr_cap <- Inf
immig_rate <- 0.05
ana_rate <- 1
sim_pars <- c(clado_rate, ext_rate, carr_cap, immig_rate, ana_rate)
set.seed(1)
island_replicates <- DAISIE_sim_cr(
  time = 1,
  M = 100,
  pars = sim_pars,
  replicates = 2,
  plot_sims = FALSE,
  verbose = FALSE
)

## Simulate 2 islands for 1 million years with two types of species (type1
## and type 2). Pool size 100
## Fraction of type 2 species in source pool is 0.15. Function will
## simulate until number of islands where type 2 species has colonised is
## equal to number specified in replicates.

clado_rate <- 0.5
ext_rate <- 0.2
carr_cap <- Inf
immig_rate <- 0.005
ana_rate <- 1
sim_pars_type1 <- c(clado_rate, ext_rate, carr_cap, immig_rate, ana_rate)
sim_pars_type2 <- sim_pars_type1 * 2
set.seed(1)
island_replicates_2types <- DAISIE_sim_cr(
  time = 1,
  M = 100,
  pars = c(sim_pars_type1, sim_pars_type2),
  replicates = 2,
  prop_type2_pool = 0.15,
  plot_sims = FALSE,
  verbose = FALSE
)
## Simulate two non-oceanic island for 1 million years.
## Pool size 500. Island area as a proportion
## of mainland is 0.1, proportion of native species is 0.9.
clado_rate <- 0.5
ext_rate <- 0.2
carr_cap <- Inf
immig_rate <- 0.005
ana_rate <- 1
sim_pars <- c(clado_rate, ext_rate, carr_cap, immig_rate, ana_rate)
set.seed(1)
island_replicates <- DAISIE_sim_cr(
  time = 1,
  M = 500,
  pars = sim_pars,
  replicates = 2,
  nonoceanic_pars = c(0.1, 0.9),
  plot_sims = FALSE,
  verbose = FALSE
)

## Simulate 2 islands for 1 million years with a shift in immigration rate
## at 0.195 Ma, and plot the species-through-time plot. Pool size 296.

pars_before_shift <- c(0.079, 0.973, Inf, 0.136, 0.413)
pars_after_shift <- c(0.079, 0.973, Inf, 0.652, 0.413)
tshift <- 0.195
set.seed(1)
island_shift_replicates <- DAISIE_sim_cr_shift(
  time = 1,
  M = 296,
  pars = c(pars_before_shift, pars_after_shift),
  replicates = 2,
  shift_times = tshift,
  plot_sims = FALSE,
  verbose = FALSE
)
```
