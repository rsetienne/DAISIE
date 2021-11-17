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
