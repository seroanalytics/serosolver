% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc_help.R
\name{setup_infection_histories_antibody_level}
\alias{setup_infection_histories_antibody_level}
\title{Propose initial infection histories based on antibody levels}
\usage{
setup_infection_histories_antibody_level(
  antibody_data,
  possible_exposure_times,
  space = 5,
  antibody_cutoff = 2,
  sample_prob = 0.9
)
}
\arguments{
\item{antibody_data}{the matrix of titres data with columns for individual, sample, and antibody level}

\item{possible_exposure_times}{vector of real times for all strains}

\item{space}{how many epochs must separate proposed infections}

\item{antibody_cutoff}{specifies how high the antibody level must be to imply an infection}

\item{sample_prob}{if antibody levels suggest an infection, then add an infection with 1 minus this probability}
}
\value{
an nxm matrix of infection histories containing 1s and 0s, where n is the number of individuals and m is the number of time periods for potential infection
}
\description{
Very similar to \code{\link{setup_infection_histories}}, but is not restricted to placing starting infections against antigens to which an individual has a measurable antibody level Given a matrix of antibody data, proposes plausible initial infection histories from which to begin MCMC sampling.
The idea is to move along time and look at an individual's antibody level against each antigen/variant Where antibodies are elevated, this suggests an infection. However, to avoid suggesting multiple infections for regions of high antigenic similarity, we place a necessary gap (defined by `space`) between proposed infection times.
}
\examples{
data(example_antibody_data)
data(example_antigenic_map)
start_inf <- setup_infection_histories_antibody_level(example_antibody_data, example_antigenic_map$inf_times)
}
\seealso{
Other setup_infection_histories: 
\code{\link{setup_infection_histories}()}
}
\concept{setup_infection_histories}
