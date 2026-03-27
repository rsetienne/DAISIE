# Create named list of trait state parameters

Create named list of trait state parameters

## Usage

``` r
create_trait_pars(
  trans_rate,
  immig_rate2,
  ext_rate2,
  ana_rate2,
  clado_rate2,
  trans_rate2,
  M2
)
```

## Arguments

- trans_rate:

  A numeric with the per capita transition rate with state1

- immig_rate2:

  A numeric with the per capita immigration rate with state2

- ext_rate2:

  A numeric with the per capita extinction rate with state2

- ana_rate2:

  A numeric with the per capita anagenesis rate with state2

- clado_rate2:

  A numeric with the per capita cladogenesis rate with state2

- trans_rate2:

  A numeric with the per capita transition rate with state2

- M2:

  A numeric with the number of species with trait state 2 on mainland

## Value

list of numerical values containing trait state parameters

## Examples

``` r
trait_pars <- create_trait_pars(
   trans_rate = 0.5,
   immig_rate2 = 0.1,
   ext_rate2 = 0.2,
   ana_rate2 = 0.3,
   clado_rate2 = 0.4,
   trans_rate2 = 0.5,
   M2 = 1000
 )
```
