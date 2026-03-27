# Computes integral of a very peaked function, modified from the SADISA package

computes the logarithm of the integral of exp(logfun) from 0 to Inf
under the following assumptions:

- "exp(logfun)" has a single, sharply peaked maximum

- "exp(logfun)" is increasing to the left of the peak and decreasing to
  the right of the peak

- "exp(logfun)" can be zero or positive at zero

- "exp(logfun)" tends to zero at infinity

## Usage

``` r
integral_peak(
  logfun,
  xx = seq(-20, 20, 2),
  xcutoff = 2,
  ymaxthreshold = 1e-12,
  pars1,
  pars2,
  brts,
  stac,
  missnumspec,
  methode,
  abstolint,
  reltolint,
  verbose,
  pick,
  par_mean,
  par_sd,
  par_upper_bound
)
```

## Arguments

- logfun:

  the logarithm of the function to integrate

- xx:

  the initial set of points on which to evaluate the function

- xcutoff:

  when the maximum has been found among the xx, this parameter sets the
  width of the interval to find the maximum in

- ymaxthreshold:

  sets the deviation allowed in finding the maximum among the xx

## Value

the result of the integration

## References

Haegeman, B. & R.S. Etienne (2017). A general sampling formula for
community structure data. Methods in Ecology & Evolution 8: 1506-1519.
https://doi.org/10.1111/2041-210X.12807
