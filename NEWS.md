# DAISIE 4.3.3

* Address problem detected by valgrind: unitialized member variable 
bulirsch_stoer<>::m_dt_last.
  * Patched version of `boost/numeric/odeint/stepper/bulirsch_stoer.hpp`. This
  is done by including the patched header file `src/patched_bulrisch_stoer.hpp` 
  this header before `boost/numeric/odeint` to shadow 
  `boost/numeriuc/odeint/stepper/bulrisch_stoer.hpp`
  The issue is *not* fixed in BOOST_VERSION == 1.81.0.
  Must check for fixes in upcomming boost (BH) releases.
* Fix tab spacing in `src/DAISIE_loglik_rhs_FORTRAN.f95`.
* Re-implement clang16 `-D_HAS_AUTO_PTR_ETC` flag fix via Makevars[.win] to 
comply with CRAN requests and to force C++ standard C++14 without using 
SystemRequirements line in DESCRIPTION, at CRAN's request.
* Skip plotting tests that cause problems in headless systems. These should
be run manually.
* Skip tests of as-of-now experimental steppers from ODEINT rosenbrock4 and 
adams bashfort moulton 1 because they are too slow.

# DAISIE 4.3.2

* Apply CRAN suggested fixes to clang16 issues with deprecated C++ functions included the Boost library, which are used in some of the stepper functions.
  * Add config.h in a macro, checking for, and setting, `_HAS_AUTO_PTR_ETC` and `BOOST_NO_AUTO_PTR`.
  * Change SystemRequirements in DESCRIPTION from C++17 to C++14.
This same fix was applied in package `'DDD'` version 5.2.1.

# DAISIE 4.3.1

* Fix issue that prevented 'covr' from running correctly.
* Depend on R >= 4.1.0, C++17.
* Tests of `DAISIE_ML_*()` and integration tests don't run on CRAN due to too long runtime.

# DAISIE 4.3.0

* Due to recent changes in CRAN policy that result in warnings due to the use of Rcpp related packages that link to libraries that use `sprintf` ('Rcpp' and 'BH'), 'DAISIE' now requires Rcpp v(>= 1.0.10) and BH v(>= 1.81.0-1). See https://github.com/RcppCore/Rcpp/pull/1236 and https://github.com/eddelbuettel/bh/pull/90 respectively. Also require package 'DDD' v(>= 5.0.0).
* Add code for the trait-dependent simulations considering two carrying capacities for each trait state.
* Improvements to the relaxed rate model fitting procedure.
* Add `DAISIE_count_species()` to count the number of species in a datalist or simulated data.
* Simply printing code by wrapping multiple instances into `print_parameters_and_loglik()`.
* Overall documentation improvements and updated references.
* Bug fix in `DAISIE_sim_cr_iw()`, which wrongly computed number of species at present. #147
* Improve the likelihood calculation. Fixed incorrect likelihood present for  type 2 case CS model inference if the values of lambda were high, in the C++ implementation.
* Overhaul GHA workflows to be more up to date and more conservative when to run. Feature branches require the tag
[run ci] be added to the commit message to run.

# DAISIE 4.2.1

* Bug fixes in `DAISIE_loglik_CS()`on the likelihood code for the continental sampling (probability of initial presence on the island).
* Bug fix in the `DAISIE_sim_cr_iw()` so replication works correctly.

# DAISIE 4.2.0

* Important bugfixes on estimation when data contains a lineage or a clade for 
which (see `stac_key` vignette for details on each case):
  * The divergence time is unknown but known to have occurred *after* a 
  specific point in time in the island's existence.
  * The divergence time is unknown but known to have occurred *before* a 
  specific point in time in the island's existence (usually known when dated 
  population level phylogenetic data is available).
* Improve and finalise the continental island estimation scenario, when there 
are initially species present on the island on inception. 
* Improve tests related to the two points above.
* Remove DAISIE time dependent estimation code for now to ease work on the remaining
code. Development will proceed separately.
* Rate computations in simulation are first calculated per lineage, to then be 
calculated per island. This will be needed for time dependent estimation later on.
* Improve zenodo metadata.
* Add `methode = "odeint::adams_bashforth_moulton_X"` to the list of available
numeric integrators for parameter estimation.


# DAISIE 4.1.1

* Correctly use `is.data.frame()` rather than `class(foo) == "data.frame"` to satisfy CRAN note.

# DAISIE 4.1.0

* Rework IW simulation code to:
    * Produces full stt tables.
    * Can run with conditioning on number of surviving species to the present.
* Added `DAISIE_CS_max_steps()` to manually change the maximum number of integration steps for the CS C++ integrators.
* Renamed and reworked simulation architecture.
* Renamed `totaltime` to `total_time`.
* Integrators from package `'deSolve'` are reinstated as the default ODE solvers.
* Added stac_key vignette detailing possible statuses of colonist lineages.
* Use LazyData.

# DAISIE 4.0.5

**N.B.: MacOS users may experience issues when installing DAISIE, especially when on MacOS Big Sur. If that is you case, please see [here](https://github.com/rsetienne/DAISIE/blob/6da0e3f65680d5f237345ef80935bda7541cf230/doc/DAISIE_macOS.md) for detailed installation instructions.**

* CS model can be fit with C++ code, similarly to the IW model since version v3.2.0. To do so, set `methode` to:
    * `"odeint::runge_kutta_cash_karp54"`
    * `"odeint::runge_kutta_fehlberg78"` [default]
    * `"odeint::runge_kutta_dopri5"`
    * `"odeint::bulirsch_stoer"`

without the `odeint::`-prefix, the `deSolve::ode` method is assumed. `deSolve::ode` using `"lsodes"` is still the default, as before. Make sure `CS_version = 0` (default) to use C++ ODEINT numerical integration.
* Better handling of number of threads to use with `DAISIE_IW_num_threads()` for ML estimation of the IW model when more threads are available.
* Fixes to numerical instability in conditional probability calculation in the IW model when the probability of extinction is very close to 1.
* Introduce an approximation when the log conditional likelihood cannot be computed and add message when this is triggered.

# DAISIE 4.0.2
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4813095.svg)](https://doi.org/10.5281/zenodo.4813095)

**N.B.: MacOS users may experience issues when installing DAISIE, especially when on MacOS Big Sur. If that is you case, please see [here](https://github.com/rsetienne/DAISIE/blob/6da0e3f65680d5f237345ef80935bda7541cf230/doc/DAISIE_macOS.md) for detailed installation instructions.**

* Suggest ggtree >= 3.0.0.

# DAISIE 4.0.1
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4813094.svg)](https://doi.org/10.5281/zenodo.4813094)

**N.B.: MacOS users may experience issues when installing DAISIE, especially when on MacOS Big Sur. If that is you case, please see [here](https://github.com/rsetienne/DAISIE/blob/6da0e3f65680d5f237345ef80935bda7541cf230/doc/DAISIE_macOS.md) for detailed installation instructions.**

* Fix possibility of fitting CS model with IW likelihood on simulated data by setting `CS_version = 0`. Improve `CS_version` documentation.

# DAISIE 4.0.0
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4813090.svg)](https://doi.org/10.5281/zenodo.4813090)

**N.B.: MacOS users may experience issues when installing DAISIE, especially when on MacOS Big Sur. If that is you case, please see [here](https://github.com/rsetienne/DAISIE/blob/6da0e3f65680d5f237345ef80935bda7541cf230/doc/DAISIE_macOS.md) for detailed installation instructions.**

* Fix bug when calculating conditional probabilities, which are now correctly calculated from island age to the present.
* Fix bug when calculating probabilities upon migration, which assumed no recolonisation was possible in the CS model. Handling recolonisation in the same manner is not possible for the IW model, so an approximation is now made. The influence of the bug and of the approximation in the IW model is expected to be minimal, particularly in cases where the colonisation rate is low. Such cases of low colonisation are the norm.
* Add `num_cycles` argument to ML functions, allowing to specify how many cycles the optimizer should take. Defaults to 1.
* Minor vignette corrections.
* Improve README.md documentation.

# DAISIE 3.2.1
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4633441.svg)](https://doi.org/10.5281/zenodo.4633441)

**N.B.: MacOS users may experience issues when installing DAISIE, especially when on MacOS Big Sur. If that is you case, please see [here](https://github.com/rsetienne/DAISIE/blob/6da0e3f65680d5f237345ef80935bda7541cf230/doc/DAISIE_macOS.md) for detailed installation instructions.**

* Minor documentation improvements.

# DAISIE 3.2.0
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4630604.svg)](https://doi.org/10.5281/zenodo.4630604)


**N.B.: MacOS users may experience issues when installing DAISIE, especially when on MacOS Big Sur. If that is you case, please see [here](https://github.com/rsetienne/DAISIE/blob/6da0e3f65680d5f237345ef80935bda7541cf230/doc/DAISIE_macOS.md) for detailed installation instructions.**

* `DAISIE_loglikg_IW()` is now more efficient and numerically stable. Numerical integration is now done via C++ with package `odeint`.
* Add relaxed rate capabilities (both inference and simulations). Relaxed rate models allow for parameters to not be static, but to be sampled by specific probability distributions.
* Introduce `MinAge` data status in DAISIE data objects. A status containing `MinAge` sets a lower boundary for colonization in situations when the precise colonization time is unknown. This is interpreted by `DAISIE_dataprep()` so that the information is passed on to the likelihood optimization functions. See the `DAISIE_dataprep()` help page for more details. In the back-end this results in new `stac` values 8 and 9.
* Bug fix of "bug 2" in the bug report manuscript. This bug was present in `DAISIE_ONEcolonist()` when recolonization occurs. It has now been fixed so that the colonization and branching times are stored in the way that we now think is the best for it to be dealt with in the likelihood code. In recolonization cases, `$other_clades_same_ancestor` renamed to `$all_colonisations`. #125
* Fix bug which occurs rarely, when computing log conditional probabilities. Only applicable to ML code running with `cond`.
* Removed deprecated legacy functions. Removed all functions named `DAISIE_*_VERSION_NUMBER()` and all `DAISIE_calc_*_rate()` funcions and `get_brts_mya()`. #126
* Made some functions internal, as they should be. `DAISIE_make_global()` and `create_island()` are now internal. #127
* @HHildenbrandt is now an author.
* Added @xieshu95's and @joshwlambert's ORCIDs.
* Added a `NEWS.md` file to track changes to the package.

# DAISIE 3.1.0
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4054059.svg)](https://doi.org/10.5281/zenodo.4054059)

* Expands the possibility of conditioning simulations and MLE of the CS model on the number of colonizing lineages. 

* In simulation and ML functions, the `cond` argument can now be greater than one. A non-zero `cond` signifies that the ML or simulation is conditioned on having at least `cond` colonizations on the island.

* Implements #121, at sim and ML level.

* Add BugReports, Website and missing ORCID in DESCRIPTION.

# DAISIE 3.0.1
* Correct @joshwlambert's name in `DESCRIPTION`.
* `DAISIE_sim_relaxed_rate()` input is closer to `DAISIE_ML()` input.
* Documentation improvements.
* Tweak Makevars.

# DAISIE 3.0.0
* Major revamp to simulation code. Simulations now accessed using `DAISIE_sim_*()` syntax.
* Constant rate, time-dependent, trait-dependent, and (multiple) split-rate available.
* Relaxed-rate inference available in `DAISIE_ML_CS()`.
* Improved vignettes documenting CS and IW cases.
* Full stt can be returned by setting `sample_freq = Inf` in `DAISIE_sim_*()` functions.
* Optional plotting with `DAISIE_plot_input()`. (Requires additional dependencies).
* Back-end architecture improvements.

# DAISIE 2.0.1
Minor update to v2.0: when empty islands are simulated the output list contains only one element instead of two (where the second indicated stac = 0, i.e. no surviving colonization).

# DAISIE 2.0

Contains the functions used in:

* Valente L., Etienne R.S., Garcia-R J.C. (2019) Deep macroevolutionary impact of humans on New Zealand's unique avifauna. Current Biology , 29, 2563-2569.e4. https://doi.org/10.1016/j.cub.2019.06.058

* Valente L., Phillimore A.B., Melo M. et al. (2020) A simple dynamic model explains the diversity of island birds worldwide.  Nature , 579, 92–96. https://doi.org/10.1038/s41586-020-2022-5

* Hauffe T, Delicado D, Etienne R.S. and Valente L. (2020) Lake expansion elevates equilibrium diversity via increasing colonisation. Journal of Biogeography. https://doi.org/10.1111/jbi.13914
