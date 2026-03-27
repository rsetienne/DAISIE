# Extract the STT median from the output of DAISIE_sim functions

Extract the STT median from the output of DAISIE_sim functions

## Usage

``` r
DAISIE_extract_stt_median(island_replicates, trait_pars = NULL)
```

## Arguments

- island_replicates:

  List output from
  [`DAISIE_sim_core_cr()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim_core_cr.md),
  [`DAISIE_sim_core_time_dep()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim_core_time_dep.md),
  [`DAISIE_sim_core_cr_shift()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim_core_cr_shift.md)
  or
  [`DAISIE_sim_min_type2()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim_min_type2.md)
  functions. Minimally, this must be a list that has as many elements as
  replicates. Each element must be a list with the elements
  `island_age`, `not_present` and `stt_all`. `stt_all` must be a data
  frame with the column names `Time`, `nI`, `nA`, `nC` and `present`.

- trait_pars:

  A named list containing diversification rates considering two trait
  states created by
  [`create_trait_pars`](https://rsetienne.github.io/DAISIE/reference/create_trait_pars.md):

  - \[1\]:A numeric with the per capita transition rate with state 1

  - \[2\]:A numeric with the per capita immigration rate with state 2

  - \[3\]:A numeric with the per capita extinction rate with state 2

  - \[4\]:A numeric with the per capita anagenesis rate with state 2

  - \[5\]:A numeric with the per capita cladogenesis rate with state 2

  - \[6\]:A numeric with the per capita transition rate with state 2

  - \[7\]:A numeric with the number of species with trait state 2 on
    mainland

## Value

a matrix (?)
