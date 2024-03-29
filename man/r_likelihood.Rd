% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/observation_model_likelihoods.R
\name{r_likelihood}
\alias{r_likelihood}
\title{Likelihood function given discrete data (normal)}
\usage{
r_likelihood(
  expected,
  data,
  theta,
  expected_indices = NULL,
  measurement_shifts = NULL
)
}
\arguments{
\item{expected}{vector of expected antibody levels}

\item{data}{vector of observed discrete measurements}

\item{expected_indices}{the indices of the measurement_shifts vector that each predicted antibody level needs adding to it}

\item{measurement_shifts}{the vector of measurement shifts for each cluster to add to the predicted antibody levels}
}
\value{
a vector with the likelihood of making each observation given the predictions
}
\description{
Calculates the likelihood of observing a set of discrete measurements given a corresponding set of predicted antibody levels
}
