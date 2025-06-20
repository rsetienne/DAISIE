% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_pars.R
\name{create_CS_version}
\alias{create_CS_version}
\title{Creates the list object for CS_version argument in DAISIE_ML_CS}
\usage{
create_CS_version(
  model = 1,
  function_to_optimize = "DAISIE",
  relaxed_par = NULL,
  par_sd = 0,
  par_upper_bound = Inf,
  integration_method = "standard",
  seed = 42,
  sample_size = 100,
  parallel = FALSE,
  n_cores = 1
)
}
\arguments{
\item{model}{the CS model to run, options are \code{1} for single rate
DAISIE model, \code{2} for multi-rate DAISIE, or \code{0} for IW test
model}

\item{function_to_optimize}{likelihood function that must be optimized in ML,
either 'DAISIE', 'DAISIE_approx', or 'DAISIE_DE'}

\item{relaxed_par}{the parameter to relax (integrate over). Options are
\code{"cladogenesis"},
\code{"extinction"},
\code{"carrying_capacity"},
\code{"immigration"},
\code{"anagenesis"}}

\item{par_sd}{standard deviation of the parameter to relax}

\item{par_upper_bound}{upper bound of the parameter to relax}

\item{integration_method}{method of integration, either 'standard','stratified'
or 'MC'}

\item{seed}{seed of the random number generator in case of 'MC'}

\item{sample_size}{size of sample in case of 'MC' or 'stratified'}

\item{parallel}{use parallel computing or not in case of 'MC' or 'stratified'}

\item{n_cores}{number of cores to use when run in parallel}
}
\value{
A list of four elements
\itemize{
  \item{model: the CS model to run, options are \code{1} for single rate
  DAISIE model, \code{2} for multi-rate DAISIE, or \code{0} for IW test
  model}
  \item{fumction_to_optimize} likelihood function that must be optimized in
  ML, either 'DAISIE', 'DAISIE_approx', or 'DAISIE_DE'
  \item{relaxed_par: the parameter to relax (integrate over), for model = 2.}
  \item{par_sd: the standard deviation of the parameter to relax}
  \item{par_upperbound: upper bound of the parameter to relax.}
  \item{integration_method: method of integration, either 'standard',
  'stratified' or 'MC'}
  \item{seed: random seed in case of integration_method = 'MC'}
  \item{sample_size: size of sample in case of integration_method = 'MC'
  or 'stratified'}
  \item{parallel: use parallel computing or not in case of integration_method
  = 'MC' or 'stratified'}
  \item{n_cores: number of cores to use when run in parallel}
}
}
\description{
Creates the list object for CS_version argument in DAISIE_ML_CS
}
