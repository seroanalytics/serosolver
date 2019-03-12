#include <Rcpp.h>
using namespace Rcpp;

#ifndef LIKELIHOOD_DATA_INDIVIDUAL_H
#define LIKELIHOOD_DATA_INDIVIDUAL_H
double likelihood_data_individual(const NumericVector &theta, 
				  const IntegerVector &infection_history, 
				  const NumericVector &circulation_times, 
				  const IntegerVector &circulation_times_indices,
				  const NumericVector &sample_times,
				  const IntegerVector &data_indices,
				  const IntegerVector &measurement_strain_indices, 
				  const NumericVector &antigenic_map_long, 
				  const NumericVector &antigenic_map_short,
				  const int &number_strains,
				  const NumericVector &data,
				  const NumericVector &to_add,
				  const double &age,
				  const Nullable<List> &additional_arguments
				  ) ;
#endif
