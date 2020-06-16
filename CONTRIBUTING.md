# DAISIE Contributing Guidelines
Feature requests and collaboration to DAISIE are encouraged. To ensure the consistency and quality of the code when so many collaborators are working on the same codebase, certain conventions must be followed. External pull-requests will be rejected or changes will be required if these conditions are not met. If this happens, please be aware that this is in no way personal, and be sure that we appreciate your interest in DAISIE, rather, we enforce these rules to ensure the package is as readable and bug-free as we can make it.

## Naming
1. Files and functions in their own file for functions start with the prefix `DAISIE_`.

2. Helper functions can go into files but do not have to start with the prefix `DAISIE_` but will be stored within a DAISIE_ prefixed file owing to point 1.

## Style
1. Code should be formatted according to the [tidyverse style guide](https://style.tidyverse.org/). This can be checked by running `devtools::lint()`, this will examine the DAISIE package to check code is in line with the style guide and give corrections if an error in style exists.

## Coverage and testing
1. Cover new functions with a unit test. Unit tests function as both a way to ensure new functions behave as expected, and as an internal form of documentation. Unit tests cover things such as: the correct type is returned when a function is called correctly, boundary cases behave as expected, and simple, short runs of the function return the expected, correct output. A helper function for unit testing is available in `usethis::use_test()`. When the cursor is on the function requiring testing, a new file with skeleton code for testing will be generated.

2. Test coverage can be calculated using the function `devtools::test_coverage()`. This will produce a window with each function’s coverage across the whole package.

3. All new functions must have (at least) 100% coverage.

4. Integration tests are encouraged. Integration tests, unlike unit tests, test the functionality of the package in a more unified way. Examples of integration testing include testing an entire ML routine. Because integration tests generally much slower to run, we request that these be skipped on local builds and only run on Travis CI. Care should still be taken to balance the test’s speed and usefulness, so do try to make it as fast as possible. If this is still impossible, a final option is to include said integration test in DAISIE’s sister package [DAISIEtesting](https://github.com/Neves-P/DAISIEtesting), built specifically for this purpose.
 
## Documentation
1. Use the R package `roxygen2` to document the package and the individual functions. Use `devtools::document()` to write the NAMESPACE. The functions are documented using the shortcut `Ctrl` + `Alt` + `Shift` + `R` while the cursor is over the function name. This will write the skeleton for the function documentation into which you can write the function description and parameters, etc.

2. Every function needs to have documentation. This includes helper or secondary functions. Internal functions that do not require an explicit documentation entry in the CRAN manual, should be marked with `@keyword Internal` to ensure that developer documentation is still accessible, while preventing unnecessary clutter for the user in the manual.

3. Parameter documentation is written in the `default_params_doc.R` and called by each function. Thus, each individual function’s roxygen block should NEVER have any `@param` field, but instead have `@inheritParams default_params_doc`.

4. All exported functions should be accompanied by runnable and easy to read examples that document all possible flows of code. This means that if a simulation can be run in the IW or CS modes, for instance, example code should be provided for both methods. In the very rare cases when running an example is not possible, e.g. it is virtually impossible to make the code fast, or an error *is* expected, please write the example anyway and encapsulate it with the `\dontthest{}` or `\dontrun{}` tags respectively.

5. Authorship/blame. Any direct or indirect contributor to a piece of code should be listed as `@author`, unless specifically requested by the author. Similarly, if code refers to a specific publication, please include a formatted citation under `@references`.

6. Exported functions. As mentioned above, the NAMESPACE is handled automatically by `roxygen2`. Only functions required by users should be exported with `@export`. For example, the `sim_core` family of functions should NOT be exported, while the `sim` family of functions is supposed to be.

## Git, GitHub and continuous integration
1. DAISIE is sorted and collaborated on via [GitHub](www.github.com), on this same [repository](www.github.com/rsetienne/DAISIE). On which you can make a branch (with the permission of the maintainer Rampal, GitHub account name: @rsetienne), this will allow you to work on code without disturbing the existing working DAISIE (`master` branch).

2. Once your code is working, it can be merged with the `develop` branch and eventually `master`. Code can be pushed to your branch using git, either through the command line or directly in RStudio. Do not merge your topic branch to `develop` if your branch fails `R CMD Check`. You can run `R CMD Check` locally using RStudio with `Ctrl` + `Shift` + `E`, on Windows you will need to have the latest version of [RTools](https://cran.r-project.org/bin/windows/Rtools/) installed to build and check the package from source. Topic branches can only be merged to `develop`, and only `develop` is merged to `master`.

3. The DAISIE package on GitHub uses continuous integration (CI) using [Travis CI](https://travis-ci.org/). This means that every time you push code to the GitHub repository it will be built, checked for errors and test coverage will be calculated and uploaded to [codecov.io](https://codecov.io/).

## Vignette
1. A vignette describes the usage and guides the user in using the package. It is useful to present the capabilities of the package to a new user. Vignettes are written in Rmarkdown and as such, comprise markdown formatted text and regular R code which is ran to generate the vignette. Vignettes are built every time the package is checked (meaning that they are rebuilt by Travis on any new pushed commit). Thus, the contained code should be readable, concise, and above all, fast to run (few replicates, high ML tolerance or iteration count, etc), otherwise `R CMD Check` will take too long and fail.
