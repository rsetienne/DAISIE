# IW concurrency control

Sets or retrieves the number of threads used by the odeint solver.

## Usage

``` r
DAISIE_IW_num_threads(num_threads)
```

## Arguments

- num_threads:

  `num_threads < 0 or omitted`: retrieves number of threads.  
  `num_threads = 0`: sets the number of threads to the number of
  available cores.  
  `num_threads = 1`: single-threaded execution.  
  `num_threads > 1`: sets the number of threads to `num_threads`.

## Value

number of threads

## Note

The maximum number of threads is limited to the value of the C++
standard library function `std::thread::hardware_concurrency()`. This is
also the default number of threads upon library load. Multithreading
incurs some overhead. Therefore, single-threaded execution might be
faster for small systems.
