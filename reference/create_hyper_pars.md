# Create list of hyperparameters

Create list of hyperparameters

## Usage

``` r
create_hyper_pars(d, x)
```

## Arguments

- d:

  Numeric defining the scaling parameter for exponent for calculating
  cladogenesis rate.

- x:

  Numeric defining the exponent for calculating extinction rate.

## Value

Named list with hyperparameters

## Author

Pedro Neves, Joshua Lambert

## Examples

``` r
hyper_pars <- create_hyper_pars(d = 0.027, x = 0.15)
```
