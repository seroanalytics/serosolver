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

#ifndef TITRE_DATA_FAST_INDIVIDUAL_BASE_H
#define TITRE_DATA_FAST_INDIVIDUAL_BASE_H
void titre_data_fast_individual_base(NumericVector &predicted_titres,
				     const double &mu, const double &mu_short, 
				     const double &wane, const double &tau,
				     const NumericVector &infection_times,
				     const IntegerVector &infection_strain_indices_tmp,
				     const IntegerVector &measurement_strain_indices,
				     const NumericVector &sample_times,
				     const int &index_in_samples,
				     const int &end_index_in_samples,
				     const int &start_index_in_data,
				     const IntegerVector &nrows_per_blood_sample,
				     const int &number_strains,
				     const NumericVector &antigenic_map_short,
				     const NumericVector &antigenic_map_long
				     );
#endif
