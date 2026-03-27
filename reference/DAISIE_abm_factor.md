# Sets or retrieves the factor to calculate the step-size used by the odeint::adams_bashforth\[\_moulton\] solvers.

Sets or retrieves the factor to calculate the step-size used by the
odeint::adams_bashforth\[\_moulton\] solvers.

## Usage

``` r
DAISIE_abm_factor(factor)
```

## Arguments

- factor:

  sets step-size to `factor * (t1 - t0)`.  

## Value

current factor
