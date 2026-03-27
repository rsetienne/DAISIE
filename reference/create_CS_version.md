# Creates the list object for CS_version argument in DAISIE_ML_CS

Creates the list object for CS_version argument in DAISIE_ML_CS

## Usage

``` r
create_CS_version(
  model = 1,
  function_to_optimize = "DAISIE",
  relaxed_par = NULL,
  par_sd = 0,
  par_upper_bound = Inf,
  integration_method = "standard",
  seed = 42,
  sample_size = 100,
  parallel = FALSE,
  n_cores = 1
)
```

## Arguments

- model:

  the CS model to run, options are `1` for single rate DAISIE model, `2`
  for multi-rate DAISIE, or `0` for IW test model

- function_to_optimize:

  likelihood function that must be optimized in ML, either 'DAISIE',
  'DAISIE_approx', or 'DAISIE_DE'

- relaxed_par:

  the parameter to relax (integrate over). Options are `"cladogenesis"`,
  `"extinction"`, `"carrying_capacity"`, `"immigration"`, `"anagenesis"`

- par_sd:

  standard deviation of the parameter to relax

- par_upper_bound:

  upper bound of the parameter to relax

- integration_method:

  method of integration, either 'standard','stratified' or 'MC'

- seed:

  seed of the random number generator in case of 'MC'

- sample_size:

  size of sample in case of 'MC' or 'stratified'

- parallel:

  use parallel computing or not in case of 'MC' or 'stratified'

- n_cores:

  number of cores to use when run in parallel

## Value

A list of four elements

- model: the CS model to run, options are `1` for single rate DAISIE
  model, `2` for multi-rate DAISIE, or `0` for IW test model

- fumction_to_optimize likelihood function that must be optimized in ML,
  either 'DAISIE', 'DAISIE_approx', or 'DAISIE_DE'

- relaxed_par: the parameter to relax (integrate over), for model = 2.

- par_sd: the standard deviation of the parameter to relax

- par_upperbound: upper bound of the parameter to relax.

- integration_method: method of integration, either 'standard',
  'stratified' or 'MC'

- seed: random seed in case of integration_method = 'MC'

- sample_size: size of sample in case of integration_method = 'MC' or
  'stratified'

- parallel: use parallel computing or not in case of integration_method
  = 'MC' or 'stratified'

- n_cores: number of cores to use when run in parallel
