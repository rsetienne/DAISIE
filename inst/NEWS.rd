\name{NEWS}
\title{NEWS}

\section{Changes in version 4.0.1}{
\strong{N.B.: MacOS users may experience issues when installing DAISIE, especially when on MacOS Big Sur. If that is you case, please see \href{https://github.com/rsetienne/DAISIE/blob/6da0e3f65680d5f237345ef80935bda7541cf230/doc/DAISIE_macOS.md}{here} for detailed installation instructions.}
\itemize{
\item Fix possibility of fitting CS model with IW likelihood on simulated data by setting \code{CS_version = 0}. Improve \code{CS_version} documentation.
}
}

\section{Changes in version 4.0.0}{
\strong{N.B.: MacOS users may experience issues when installing DAISIE, especially when on MacOS Big Sur. If that is you case, please see \href{https://github.com/rsetienne/DAISIE/blob/6da0e3f65680d5f237345ef80935bda7541cf230/doc/DAISIE_macOS.md}{here} for detailed installation instructions.}
\itemize{
\item Fix bug when calculating conditional probabilities, which are now correctly calculated from island age to the present.
\item Fix bug when calculating probabilities upon migration, which assumed no recolonisation was possible in the CS model. Handling recolonisation in the same manner is not possible for the IW model, so an approximation is now made. The influence of the bug and of the approximation in the IW model is expected to be minimal, particularly in cases where the colonisation rate is low. Such cases of low colonisation are the norm.
\item Add \code{num_cycles} argument to ML functions, allowing to specify how many cycles the optimizer should take. Defaults to 1.
\item Minor vignette corrections.
\item Improve README.md documentation.
}
}

\section{Changes in version 3.2.1}{
\href{https://doi.org/10.5281/zenodo.4633441}{\figure{https://zenodo.org/badge/DOI/10.5281/zenodo.4633441.svg}}

\strong{N.B.: MacOS users may experience issues when installing DAISIE, especially when on MacOS Big Sur. If that is you case, please see \href{https://github.com/rsetienne/DAISIE/blob/6da0e3f65680d5f237345ef80935bda7541cf230/doc/DAISIE_macOS.md}{here} for detailed installation instructions.}
\itemize{
\item Minor documentation improvements.
}
}

\section{Changes in version 3.2.0}{
\href{https://doi.org/10.5281/zenodo.4630604}{\figure{https://zenodo.org/badge/DOI/10.5281/zenodo.4630604.svg}}

\strong{N.B.: MacOS users may experience issues when installing DAISIE, especially when on MacOS Big Sur. If that is you case, please see \href{https://github.com/rsetienne/DAISIE/blob/6da0e3f65680d5f237345ef80935bda7541cf230/doc/DAISIE_macOS.md}{here} for detailed installation instructions.}
\itemize{
\item \code{DAISIE_loglikg_IW()} is now more efficient and numerically stable. Numerical integration is now done via C++ with package \code{odeint}.
\item Add relaxed rate capabilities (both inference and simulations). Relaxed rate models allow for parameters to not be static, but to be sampled by specific probability distributions.
\item Introduce \code{MinAge} data status in DAISIE data objects. A status containing \code{MinAge} sets a lower boundary for colonization in situations when the precise colonization time is unknown. This is interpreted by \code{DAISIE_dataprep()} so that the information is passed on to the likelihood optimization functions. See the \code{DAISIE_dataprep()} help page for more details. In the back-end this results in new \code{stac} values 8 and 9.
\item Bug fix of "bug 2" in the bug report manuscript. This bug was present in \code{DAISIE_ONEcolonist()} when recolonization occurs. It has now been fixed so that the colonization and branching times are stored in the way that we now think is the best for it to be dealt with in the likelihood code. In recolonization cases, \verb{$other_clades_same_ancestor} renamed to \verb{$all_colonisations}. #125
\item Fix bug which occurs rarely, when computing log conditional probabilities. Only applicable to ML code running with \code{cond}.
\item Removed deprecated legacy functions. Removed all functions named \verb{DAISIE_*_VERSION_NUMBER()} and all \verb{DAISIE_calc_*_rate()} funcions and \code{get_brts_mya()}. #126
\item Made some functions internal, as they should be. \code{DAISIE_make_global()} and \code{create_island()} are now internal. #127
\item @HHildenbrandt is now an author.
\item Added @xieshu95's and @joshwlambert's ORCIDs.
\item Added a \code{NEWS.md} file to track changes to the package.
}
}

\section{Changes in version 3.1.0}{
\href{https://doi.org/10.5281/zenodo.4054059}{\figure{https://zenodo.org/badge/DOI/10.5281/zenodo.4054059.svg}}
\itemize{
\item Expands the possibility of conditioning simulations and MLE of the CS model on the number of colonizing lineages.
\item In simulation and ML functions, the \code{cond} argument can now be greater than one. A non-zero \code{cond} signifies that the ML or simulation is conditioned on having at least \code{cond} colonizations on the island.
\item Implements #121, at sim and ML level.
\item Add BugReports, Website and missing ORCID in DESCRIPTION.
}
}

\section{Changes in version 3.0.1}{
\itemize{
\item Correct @joshwlambert's name in \code{DESCRIPTION}.
\item \code{DAISIE_sim_relaxed_rate()} input is closer to \code{DAISIE_ML()} input.
\item Documentation improvements.
\item Tweak Makevars.
}
}

\section{Changes in version 3.0.0}{
\itemize{
\item Major revamp to simulation code. Simulations now accessed using \verb{DAISIE_sim_*()} syntax.
\item Constant rate, time-dependent, trait-dependent, and (multiple) split-rate available.
\item Relaxed-rate inference available in \code{DAISIE_ML_CS()}.
\item Improved vignettes documenting CS and IW cases.
\item Full stt can be returned by setting \code{sample_freq = Inf} in \verb{DAISIE_sim_*()} functions.
\item Optional plotting with \code{DAISIE_plot_input()}. (Requires additional dependencies).
\item Back-end architecture improvements.
}
}

\section{Changes in version 2.0.1}{
Minor update to v2.0: when empty islands are simulated the output list contains only one element instead of two (where the second indicated stac = 0, i.e. no surviving colonization).
}

\section{Changes in version 2.0}{
Contains the functions used in:
\itemize{
\item Valente L., Etienne R.S., Garcia-R J.C. (2019) Deep macroevolutionary impact of humans on New Zealand's unique avifauna. Current Biology , 29, 2563-2569.e4. https://doi.org/10.1016/j.cub.2019.06.058
\item Valente L., Phillimore A.B., Melo M. et al. (2020) A simple dynamic model explains the diversity of island birds worldwide.  Nature , 579, 92–96. https://doi.org/10.1038/s41586-020-2022-5
\item Hauffe T, Delicado D, Etienne R.S. and Valente L. (2020) Lake expansion elevates equilibrium diversity via increasing colonisation. Journal of Biogeography. https://doi.org/10.1111/jbi.13914
}
}

