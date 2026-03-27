# Samples what event to happen next

Samples what event to happen next

## Usage

``` r
DAISIE_sample_event_cr(rates)
```

## Arguments

- rates:

  named list of numeric rates as returned by
  [`update_rates()`](https://rsetienne.github.io/DAISIE/reference/update_rates.md).

## Value

numeric indicating what event will happen, or a supposed event that
would happen in some timesteps of the ontogeny algorithm.

- \[1\]: immigration event

- \[2\]: extinction event

- \[3\]: cladogenesis event

- \[4\]: anagenesis event

## Author

Pedro Neves
