# DAISIE: Dynamic Assembly of Island biota through Speciation, Immigration and Extinction

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/DAISIE)](https://cran.r-project.org/package=DAISIE)
[![](http://cranlogs.r-pkg.org/badges/grand-total/DAISIE)]( https://CRAN.R-project.org/package=DAISIE)
[![](http://cranlogs.r-pkg.org/badges/DAISIE)](https://CRAN.R-project.org/package=DAISIE)

Branch|[Travis](https://travis-ci.org)|[Codecov](https://www.codecov.io)
---|---|---
`master`|[![Build Status](https://travis-ci.org/rsetienne/DAISIE.svg?branch=master)](https://travis-ci.org/rsetienne/DAISIE)|[![codecov.io](https://codecov.io/github/rsetienne/DAISIE/coverage.svg?branch=master)](https://codecov.io/github/rsetienne/DAISIE/branch/master)
`develop`|[![Build Status](https://travis-ci.org/rsetienne/DAISIE.svg?branch=develop)](https://travis-ci.org/rsetienne/DAISIE)|[![codecov.io](https://codecov.io/github/rsetienne/DAISIE/coverage.svg?branch=develop)](https://codecov.io/github/rsetienne/DAISIE/branch/develop)
`rampal`|[![Build Status](https://travis-ci.org/rsetienne/DAISIE.svg?branch=rampal)](https://travis-ci.org/rsetienne/DAISIE)|[![codecov.io](https://codecov.io/github/rsetienne/DAISIE/coverage.svg?branch=rampal)](https://codecov.io/github/rsetienne/DAISIE/branch/rampal)
`luis`|[![Build Status](https://travis-ci.org/rsetienne/DAISIE.svg?branch=luis)](https://travis-ci.org/rsetienne/DAISIE)|[![codecov.io](https://codecov.io/github/rsetienne/DAISIE/coverage.svg?branch=luis)](https://codecov.io/github/rsetienne/DAISIE/branch/luis)


## Installing DAISIE

The DAISIE package has a stable version on [CRAN](https://CRAN.R-project.org/package=DAISIE) and
a development version on GitHub.

### From CRAN

From within R, do:

```
install.packages("DAISIE")
```

### From GitHub

Install DAISIE from this GitHub repository by running:

```
install.packages("remotes")
remotes::install_github("rsetienne/DAISIE")
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

 * `master`: build should always pass. [@rsetienne](https://github.com/rsetienne) has control over `develop` to `master` merges.
 * `develop`: merge of topic branches, merge with `master` by [@rsetienne](https://github.com/rsetienne) iff build passes.
 * `rampal`: [@rsetienne's](https://github.com/rsetienne) personal branch.
 * `luis`: [@luislvalente's](https://github.com/luislvalente) personal branch.

## Contributors

DAISIE was originally developed by Rampal S. Etienne, Luis Valente, Albert B. Phillimore and Bart Haegeman.

Additional members working on expanding DAISIE at the [TECE lab](https://github.com/tece-lab), University of Groningen are:
Joshua Lambert, Pedro S. Neves, Richel J. C. Bilderbeek, Sebastian Mader, Shu Xie.

## References

Etienne R. S., Valente, L., Phillimore, A. B., Haegeman, B., Lambert, J. W., Neves, P., Xie, S., & Bilderbeek, R. J. C. (2020). DAISIE: Dynamical Assembly of Islands by Speciation, Immigration and Extinction. R package version 3.0.0. https://cran.r-project.org/package=DAISIE

Valente, L.M., Etienne, R.S., & Phillimore, A.B. (2014). The effects of island ontogeny on species diversity and phylogeny. Proceedings of the Royal Society B: Biological Sciences, 281(1784), 20133227–20133227. http://doi.org/10.1098/rspb.2013.3227

Valente, L.M., Phillimore, A.B., & Etienne, R.S. (2015). Equilibrium and non-equilibrium dynamics simultaneously operate in the Galápagos islands. Ecology Letters, 18(8), 844–852. http://doi.org/10.1111/ele.12461

Valente, L.M., Etienne, R.S., & Dávalos, L.M. (2017). Recent extinctions disturb path to equilibrium diversity in Caribbean bats. Nature Ecology & Evolution, 1(2), 0026. http://doi.org/10.1038/s41559-016-0026

Valente, L., Illera, J.C., Havenstein, K., Pallien, T., Etienne, R.S., & Tiedemann, R. (2017). Equilibrium Bird Species Diversity in Atlantic Islands. Current Biology, 27(11), 1660-1666. https://doi.org/(...)16/j.cub.2017.04.053

Valente, L., Phillimore, A.B., & Etienne, R.S. (2018). Using molecular phylogenies in island biogeography: It’s about time. Ecography, 1–3. http://doi.org/10.1111/ecog.03503

Valente, L., Etienne, R.S., & Garcia-R., J.C. (2019). “Deep Macroevolutionary Impact of Humans on New Zealand’s Unique Avifauna.” Current Biology 29 (15): 2563-2569.e4. https://doi.org/10.1016/j.cub.2019.06.058.

Valente, L., Phillimore, A.B., Melo, M., Warren, B.H., Clegg, S.M., Havenstein, K., Tiedemann, R., Illera, J.C.,, Thebaud, C., Aschenbach, T. & Etienne, R.S. (2020). “A Simple Dynamic Model Explains the Diversity of Island Birds Worldwide.” Nature 579 (7797): 92–96. https://doi.org/10.1038/s41586-020-2022-5.

Hauffe, T., Delicado, D., Etienne, R.S., & Valente, L. (2020). Lake expansion elevates equilibrium diversity via increasing colonization. Journal of Biogeography. In press. https://doi.org/10.1111/jbi.13914.
