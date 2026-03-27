# Samples what event to happen next

Samples what event to happen next

## Usage

``` r
DAISIE_sample_event_trait_dep(rates)
```

## Arguments

- rates:

  named list of numeric rates as returned by
  [`update_rates()`](https://rsetienne.github.io/DAISIE/reference/update_rates.md).

## Value

numeric indicating what event will happen, or a supposed event that
would happen in some timesteps of the ontogeny algorithm.

- \[1\]: immigration event with trait1

- \[2\]: extinction event with trait1

- \[3\]: cladogenesis event with trait1

- \[4\]: anagenesis event with trait1

- \[5\]: transition event with trait1

- \[6\]: immigration event with trait2

- \[7\]: extinction event with trait2

- \[8\]: cladogenesis event with trait2

- \[9\]: anagenesis event with trait2

- \[10\]: transition event with trait2

## Author

Shu Xie
