#ifndef WANE_FUNCTION_H
#define WANE_FUNCTION_H

#include <Rcpp.h>
using namespace Rcpp;

// Function to calculate non-linear waning
inline double waning_function(NumericVector theta, double time_infected, double wane);

#endif