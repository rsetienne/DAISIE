# Create vector of model parameters

Create vector of model parameters

## Usage

``` r
create_pars(clado_rate, ext_rate, carr_cap, immig_rate, ana_rate)
```

## Arguments

- clado_rate:

  Numeric rate of cladogenesis

- ext_rate:

  Numeric rate of extinction

- carr_cap:

  Numeric carrying capacity

- immig_rate:

  Numeric rate of immigration

- ana_rate:

  Numeric rate of anagenesis

## Value

Numeric vector with 5 elements

## Examples

``` r
create_pars(clado_rate = 1,
            ext_rate = 1,
            carr_cap = 10,
            immig_rate = 0.1,
            ana_rate = 1)
#> [1]  1.0  1.0 10.0  0.1  1.0
```
