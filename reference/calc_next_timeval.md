# Calculates when the next timestep will be.

Calculates when the next timestep will be.

## Usage

``` r
calc_next_timeval(max_rates, timeval, total_time)
```

## Arguments

- max_rates:

  named list of max rates as returned by
  [`update_rates`](https://rsetienne.github.io/DAISIE/reference/update_rates.md).

- timeval:

  current time of simulation

## Value

named list with numeric vector containing the time of the next timestep
and the change in time.

## Author

Joshua Lambert, Pedro Neves, Shu Xie
