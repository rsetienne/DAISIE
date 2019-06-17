# DAISIE

Branch|[Travis](https://travis-ci.org)|[Codecov](https://www.codecov.io)
---|---|---
`master`|[![Build Status](https://travis-ci.org/rsetienne/DAISIE.svg?branch=master)](https://travis-ci.org/rsetienne/DAISIE)|[![codecov.io](https://codecov.io/github/rsetienne/DAISIE/coverage.svg?branch=master)](https://codecov.io/github/rsetienne/DAISIE/branch/master)
`develop`|[![Build Status](https://travis-ci.org/rsetienne/DAISIE.svg?branch=develop)](https://travis-ci.org/rsetienne/DAISIE)|[![codecov.io](https://codecov.io/github/rsetienne/DAISIE/coverage.svg?branch=develop)](https://codecov.io/github/rsetienne/DAISIE/branch/develop)
`richel_mainland_extinction`|[![Build Status](https://travis-ci.org/rsetienne/DAISIE.svg?branch=richel_mainland_extinction)](https://travis-ci.org/rsetienne/DAISIE)|[![codecov.io](https://codecov.io/github/rsetienne/DAISIE/coverage.svg?branch=richel_mainland_extinction)](https://codecov.io/github/rsetienne/DAISIE/branch/richel_mainland_extinction)
`pedro_ontogeny`|[![Build Status](https://travis-ci.org/rsetienne/DAISIE.svg?branch=pedro_ontogeny)](https://travis-ci.org/rsetienne/DAISIE)|[![codecov.io](https://codecov.io/github/rsetienne/DAISIE/coverage.svg?branch=pedro_ontogeny)](https://codecov.io/github/rsetienne/DAISIE/branch/pedro_ontogeny)
`sebastian_archipelago`|[![Build Status](https://travis-ci.org/rsetienne/DAISIE.svg?branch=sebastian_archipelago)](https://travis-ci.org/rsetienne/DAISIE)|[![codecov.io](https://codecov.io/github/rsetienne/DAISIE/coverage.svg?branch=sebastian_archipelago)](https://codecov.io/github/rsetienne/DAISIE/branch/sebastian_archipelago)
`shu_traits`|[![Build Status](https://travis-ci.org/rsetienne/DAISIE.svg?branch=shu_traits)](https://travis-ci.org/rsetienne/DAISIE)|[![codecov.io](https://codecov.io/github/rsetienne/DAISIE/coverage.svg?branch=shu_traits)](https://codecov.io/github/rsetienne/DAISIE/branch/shu_traits)

Dynamic Assembly of Island biota through Speciation, Immigration and Extinction in R

This is a development version before the official release on CRAN.

## Installing DAISIE

The DAISIE package has a stable version on CRAN and
a development version on GitHub.

### From CRAN

From within R, do:

```
install.packages("DAISIE")
```

### From GitHub

Because the DAISIE package is located in the folder `DAISIE`, do:

```
devtools::install_github("rsetienne/DAISIE")
```

## Using DAISIE as a package dependency

### From CRAN

To your DESCRIPTION file, add `DAISIE` as any normal package.

If your package directly uses `DAISIE`:

```
Imports:
  DAISIE
```

If your package uses `DAISIE` in its peripherals (e.g. vignettes and tests):

```
Suggests:
  DAISIE
```

### From GitHub

```
Remotes:
  rsetienne/DAISIE
```

## `git` branching model

 * `master`: build should always pass. [@rsetienne](htpps://github.com/rsetienne) has control over `develop` to `master` merges.
 * `develop`: merge of topic branches, merge with `master` by [@rsetienne](htpps://github.com/rsetienne) iff build passes.
 * `pedro_ontogeny`: [@Neves-P's](https://github.com/Neves-P) topic branch adding island ontongeny functionality for simulation and parameter estimation.
 * `richel_mainland_extinction`: [@richelbilderbeek's](https://github.com/richelbilderbeek) topic branch adding dynamic mainland processes.
 * `shu_traits`: [@xieshu95's](https://github.com/xieshu95) topic branch adding lineage trait functionality
 * `sebastian_archipelago`: [@sebmader](https://github.com/sebmader) topic branch adding archipelago dynamics to DAISIE in C++
 
[@Neves-P's](https://github.com/Neves-P) is responsible for day-to-day merging and repository maintenance. For questions open an issue or contact [@Neves-P's](https://github.com/Neves-P).

## Contributors

DAISIE was originally developed by Rampal S. Etienne, Luis Valente, Albert B. Phillimore and Bart Haegeman.

Additional members working on expanding DAISIE at the [TECE lab](https://github.com/tece-lab), University of Groningen are:
Joshua Lambert, Pedro S. Neves, Richel J. C. Bilderbeek, Sebastian Mader, Shu Xie

## References

Etienne, R. S., Valente, L. M., Phillimore, A. B., & Haegeman, B., (2018). DAISIE: Dynamical Assembly of Islands by Speciation, Immigration and Extinction. R package version 1.6.0.
  
Valente, L. M., Etienne, R. S., & Phillimore, A. B. (2014). The effects of island ontogeny on species diversity and phylogeny. Proceedings of the Royal Society B: Biological Sciences, 281(1784), 20133227–20133227. http://doi.org/10.1098/rspb.2013.3227

Valente, L. M., Phillimore, A. B., & Etienne, R. S. (2015). Equilibrium and non-equilibrium dynamics simultaneously operate in the Galápagos islands. Ecology Letters, 18(8), 844–852. http://doi.org/10.1111/ele.12461

Valente, L. M., Etienne, R. S., & Dávalos, L. M. (2017). Recent extinctions disturb path to equilibrium diversity in Caribbean bats. Nature Ecology & Evolution, 1(2), 0026. http://doi.org/10.1038/s41559-016-0026

Valente, L., Phillimore, A. B., & Etienne, R. S. (2018). Using molecular phylogenies in island biogeography: It’s about time. Ecography, 1–3. http://doi.org/10.1111/ecog.03503

