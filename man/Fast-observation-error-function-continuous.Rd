% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{Fast observation error function continuous}
\alias{Fast observation error function continuous}
\alias{likelihood_func_fast_continuous}
\alias{likelihood_func_fast_continuous_fp}
\title{Fast observation error function continuous
 Calculate the probability of a set of observed antibody levels given a corresponding set of predicted antibody levels assuming continuous, bounded observations.}
\usage{
likelihood_func_fast_continuous(theta, obs, predicted_antibody_levels)

likelihood_func_fast_continuous_fp(theta, obs, predicted_antibody_levels)
}
\arguments{
\item{theta}{NumericVector, a named parameter vector giving the normal distribution standard deviation and the max observable antibody level. 
Also a parameter fp_rate, giving the probability of a (uniformly distributed) false positive given true negative.}

\item{obs}{NumericVector, the vector of observed log antibody levels}

\item{predicted_antibody_levels}{NumericVector, the vector of predicted log antibody levels}

\item{a}{vector of same length as the input data giving the probability of observing each observation given the predictions}
}
\value{
a likelihood for each observed antibody level

a likelihood for each observed antibody level
}
\description{
Fast observation error function continuous
 Calculate the probability of a set of observed antibody levels given a corresponding set of predicted antibody levels assuming continuous, bounded observations.

Fast observation error function continuous with false positives
 Calculate the probability of a set of observed antibody levels given a corresponding set of predicted antibody levels assuming continuous, bounded observations. For true negatives (i.e., model predicts no infections), then the majority of the PDF is at min_measurement. There is a probability, fp_rate, of observing a value within the detectable range.
}
\seealso{
Other likelihood_functions: 
\code{\link{likelihood_func_fast}()}

Other likelihood_functions: 
\code{\link{likelihood_func_fast}()}
}
\concept{likelihood_functions}
