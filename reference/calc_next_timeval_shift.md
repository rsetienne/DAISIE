# Calculates when the next timestep will be, and if a shift has occured.

Calculates when the next timestep will be, and if a shift has occured.

## Usage

``` r
calc_next_timeval_shift(max_rates, timeval, dynamic_shift_times, total_time)
```

## Arguments

- max_rates:

  named list of max rates as returned by
  [`update_rates`](https://rsetienne.github.io/DAISIE/reference/update_rates.md).

- timeval:

  current time of simulation

- dynamic_shift_times:

  numeric vector of times of rate shifts.

## Value

named list with numeric vector containing the time of the next timestep
and the change in time.

## Author

Joshua Lambert, Pedro Neves, Shu Xie
