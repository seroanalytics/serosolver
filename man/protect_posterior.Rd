% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{protect_posterior}
\alias{protect_posterior}
\title{Protect function (posterior function)}
\usage{
protect_posterior(f)
}
\arguments{
\item{f}{the function to be protected}
}
\value{
the protected function
}
\description{
Wrapper function to protect calls to the posterior function. If posterior does not compute correctly, returns -100000.
}
